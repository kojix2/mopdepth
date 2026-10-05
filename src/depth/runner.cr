require "hts"
require "./config"
require "./core/coverage_calculator"
require "./core/cigar"
require "./core/coverage"
require "./core/target"
require "./io/bed_reader"
require "./io/output_manager"
require "./stats/depth_stat"
require "./stats/int_histogram"
require "./stats/distribution"
require "./stats/quantize"

module Depth
  class Runner
    extend Core::CoverageUtils
    extend Stats::Distribution

    def initialize(@config : Config)
    end

    def run
      bam = nil.as(HTS::Bam?)
      output = nil.as(FileIO::OutputManager?)
      primary_error = nil.as(Exception?)
      completed = false
      begin
        @config.validate!
        opts = @config.to_options
        bam = open_bam
        all_targets = load_targets(bam)
        region = FileIO.parse_region_str(@config.chrom, all_targets)
        targets = selected_targets(all_targets, region)
        bed_map = load_bed_map(all_targets)
        output = FileIO::OutputManager.new(@config)
        write_threshold_header(output)
        state = ProcessingState.new(@config.use_median?)
        process_targets(bam, opts, targets, region, bed_map, output, state)
        write_total_outputs(output, state)
        completed = true
      rescue ex
        primary_error = ex
      ensure
        cleanup_errors = [] of String
        if manager = output
          begin
            manager.close_all(build_indices: completed)
          rescue ex
            cleanup_errors << ex.message.to_s
          end
        end
        if input = bam
          begin
            input.close
          rescue ex
            cleanup_errors << "failed to close input: #{ex.message}"
          end
        end

        if original = primary_error
          if cleanup_errors.empty?
            raise original
          else
            raise RunError.new("#{original.message}; cleanup failures: #{cleanup_errors.join("; ")}")
          end
        elsif !cleanup_errors.empty?
          raise OutputError.new(cleanup_errors.join("; "))
        end
      end
    end

    private class ProcessingState
      getter global_dist : Array(Int64)
      getter total_global_dist : Array(Int64)
      getter region_dist : Array(Int64)
      getter total_region_dist : Array(Int64)
      getter cs : Stats::IntHistogram
      property global_stat : Stats::DepthStat
      property global_region_stat : Stats::DepthStat

      def initialize(use_median : Bool)
        @global_dist = Array(Int64).new(512, 0_i64)
        @total_global_dist = Array(Int64).new(512, 0_i64)
        @region_dist = Array(Int64).new(512, 0_i64)
        @total_region_dist = Array(Int64).new(512, 0_i64)
        @global_stat = Stats::DepthStat.new
        @global_region_stat = Stats::DepthStat.new
        @cs = Stats::IntHistogram.new(use_median ? 65_536 : 0)
      end
    end

    private struct TargetSlice
      getter query_region : Core::Region
      getter offset : Int32
      getter effective_len : Int32
      getter target_size : Int32

      def initialize(@query_region : Core::Region, @offset : Int32, @effective_len : Int32)
        @target_size = effective_len + 1
      end
    end

    private def open_bam : HTS::Bam
      bam = HTS::Bam.open(@config.path, fai: @config.fasta, threads: @config.threads)
      begin
        bam.load_index
      rescue ex
        bam.close rescue nil
        raise "Failed to load index for #{@config.path}: #{ex.message}"
      end
      bam
    end

    private def write_threshold_header(output : FileIO::OutputManager)
      return unless @config.has_thresholds?

      output.write_thresholds_header(@config.threshold_values)
    end

    private def load_targets(bam : HTS::Bam) : Array(Core::Target)
      target_names = bam.header.target_names
      target_lengths = bam.header.target_len
      target_names.map_with_index do |name, i|
        length = target_lengths[i].to_i64
        unless length > 0 && length < Int32::MAX
          raise ConfigError.new("Reference length is outside the supported range for #{name}: #{length}")
        end
        Core::Target.new(name, length.to_i32, i)
      end
    end

    private def selected_targets(targets : Array(Core::Target), region : Core::Region?) : Array(Core::Target)
      return targets unless selected_region = region

      selected = targets.select { |target| target.name == selected_region.chrom }
      raise ConfigError.new("Chromosome not found: #{selected_region.chrom}") if selected.empty?
      selected
    end

    private def load_bed_map(targets : Array(Core::Target)) : Hash(String, Array(Core::Region))?
      return unless bed_path = @config.bed_path

      FileIO.read_bed(bed_path, targets)
    end

    private def process_targets(bam : HTS::Bam, opts : Core::Options, targets : Array(Core::Target),
                                region : Core::Region?, bed_map : Hash(String, Array(Core::Region))?,
                                output : FileIO::OutputManager, state : ProcessingState)
      calculator = Core::CoverageCalculator.new(bam, opts)
      coverage = Core::Coverage.new(0)
      window = @config.window_size

      targets.each do |target|
        next if skip_target?(target, bed_map)

        slice = target_slice(target, region, opts.fragment_mode)
        prepare_coverage!(coverage, slice.target_size)
        tid = calculator.calculate(coverage, slice.query_region, slice.offset)
        next if tid == Core::CoverageResult::ChromNotFound.value

        finalize_coverage!(coverage, tid, slice.target_size)
        process_target_outputs(target, coverage, tid, slice, window, bed_map, output, state)
      end
    end

    private def skip_target?(target : Core::Target, bed_map : Hash(String, Array(Core::Region))?) : Bool
      return false unless @config.no_per_base? && bed_map
      return false if @config.has_quantize? || @config.has_thresholds?

      !bed_map.has_key?(target.name)
    end

    private def target_slice(target : Core::Target, region : Core::Region?, fragment_mode : Bool) : TargetSlice
      full_region = Core::Region.new(target.name, 0, 0)
      return TargetSlice.new(full_region, 0, target.length) unless selected_region = region
      return TargetSlice.new(full_region, 0, target.length) if selected_region.start == 0 && selected_region.stop == 0

      offset = selected_region.start
      effective_len = selected_region.stop - offset
      query_region = fragment_mode ? full_region : selected_region

      # A fragment can overlap the selected interval even when neither read does.
      # Query the whole chromosome and clip its fragment events into the local buffer.
      TargetSlice.new(query_region, offset, effective_len)
    end

    private def prepare_coverage!(coverage : Core::Coverage, target_size : Int32)
      if coverage.size < target_size
        coverage.concat(Array(Int32).new(target_size - coverage.size, 0))
      end
      clear_coverage!(coverage, target_size)
    end

    private def clear_coverage!(coverage : Core::Coverage, target_size : Int32)
      i_full = 0
      while i_full < target_size
        coverage[i_full] = 0
        i_full += 1
      end
    end

    private def finalize_coverage!(coverage : Core::Coverage, tid : Int32, target_size : Int32)
      return if tid == Core::CoverageResult::NoData.value

      i = 0
      sum = 0
      while i < target_size
        sum += coverage[i]
        coverage[i] = sum
        i += 1
      end
    end

    private def process_target_outputs(target : Core::Target, coverage : Core::Coverage, tid : Int32,
                                       slice : TargetSlice,
                                       window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                       output : FileIO::OutputManager, state : ProcessingState)
      write_per_base_intervals(target, coverage, tid, slice, output)
      write_quantized_intervals(target, coverage, tid, output, slice.offset, slice.target_size) if output.f_quantized && @config.has_quantize?

      chrom_region_stat = write_regions_if_needed(target, coverage, tid, slice, window, bed_map, output, state)
      write_target_stats(target, coverage, tid, slice.target_size, chrom_region_stat, output, state)
      write_target_distributions(target, output, state)
    end

    private def write_per_base_intervals(target : Core::Target, coverage : Core::Coverage, tid : Int32,
                                         slice : TargetSlice, output : FileIO::OutputManager)
      return unless output.f_perbase

      if tid == Core::CoverageResult::NoData.value
        output.write_per_base_interval(target.name, slice.offset, slice.offset + slice.effective_len, 0)
      else
        self.class.each_constant_segment(coverage, slice.target_size - 1) do |(s, e, v)|
          output.write_per_base_interval(target.name, s + slice.offset, e + slice.offset, v)
        end
      end
    end

    private def write_regions_if_needed(target : Core::Target, coverage : Core::Coverage, tid : Int32,
                                        slice : TargetSlice,
                                        window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                        output : FileIO::OutputManager, state : ProcessingState) : Stats::DepthStat
      return Stats::DepthStat.new unless output.f_regions

      write_region_stats(target, coverage, tid, window, bed_map, state.cs, output,
        state.region_dist, slice.offset, slice.effective_len)
    end

    private def write_target_stats(target : Core::Target, coverage : Core::Coverage, tid : Int32, target_size : Int32,
                                   chrom_region_stat : Stats::DepthStat, output : FileIO::OutputManager,
                                   state : ProcessingState)
      return if tid == Core::CoverageResult::NoData.value

      self.class.bump_distribution!(state.global_dist, coverage, 0, target_size - 1)
      chrom_stat = Stats::DepthStat.from_array(coverage, 0, target_size - 2)
      state.global_stat = state.global_stat + chrom_stat
      output.write_summary_line(target.name, chrom_stat)
      write_target_region_stat(target, chrom_region_stat, output, state)
    end

    private def write_target_region_stat(target : Core::Target, chrom_region_stat : Stats::DepthStat,
                                         output : FileIO::OutputManager, state : ProcessingState)
      return unless output.f_regions

      state.global_region_stat = state.global_region_stat + chrom_region_stat
      output.write_summary_line("#{target.name}_region", chrom_region_stat)
    end

    private def write_target_distributions(target : Core::Target, output : FileIO::OutputManager, state : ProcessingState)
      if f_global = output.f_global
        self.class.write_distribution(f_global.as(::IO), target.name, state.global_dist)
        self.class.sum_into!(state.total_global_dist, state.global_dist)
      end
      if f_region = output.f_region
        self.class.write_distribution(f_region.as(::IO), target.name, state.region_dist)
        self.class.sum_into!(state.total_region_dist, state.region_dist)
        state.region_dist.fill(0_i64)
      end
      state.global_dist.fill(0_i64)
    end

    private def write_total_outputs(output : FileIO::OutputManager, state : ProcessingState)
      output.write_summary_total(state.global_stat)
      output.write_summary_line("total_region", state.global_region_stat) if output.f_regions
      if f_global = output.f_global
        self.class.write_distribution(f_global.as(::IO), "total", state.total_global_dist)
      end
      if f_region = output.f_region
        self.class.write_distribution(f_region.as(::IO), "total", state.total_region_dist)
      end
    end

    private def write_region_depth_value(output : FileIO::OutputManager, chrom : String,
                                         start : Int32, stop : Int32, name : String?,
                                         tid : Int32, value : Float64, sum : UInt64, length : Int32)
      if tid == Core::CoverageResult::NoData.value
        output.write_region_zero(chrom, start, stop, name)
      elsif @config.use_median?
        output.write_region_stat(chrom, start, stop, name, value)
      else
        output.write_region_mean(chrom, start, stop, name, sum, length)
      end
    end

    private def write_region_stats(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                   window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                   cs : Stats::IntHistogram, output : FileIO::OutputManager,
                                   region_dist : Array(Int64), offset : Int32, effective_len : Int32) : Stats::DepthStat
      if window > 0
        process_window_regions(t, coverage, tid, window, cs, output, region_dist, offset, effective_len)
      else
        process_bed_regions(t, coverage, tid, bed_map, cs, output, region_dist, offset, effective_len)
      end
    end

    private def process_window_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                       window : Int32, cs : Stats::IntHistogram,
                                       output : FileIO::OutputManager, region_dist : Array(Int64),
                                       offset : Int32, effective_len : Int32) : Stats::DepthStat
      chrom_region_stat = Stats::DepthStat.new
      start_local = 0
      end_local = effective_len
      while start_local < end_local
        stop_local = Math.min(start_local + window, end_local)
        start_abs = offset + start_local
        stop_abs = offset + stop_local
        me = 0.0
        mean_sum = 0_u64
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median?
            cs.clear
            (start_local...stop_local).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            len = stop_local - start_local
            (start_local...stop_local).each do |i|
              depth = coverage[i]
              mean_sum += depth.to_u64 if depth > 0
            end
            me = len > 0 ? mean_sum.to_f / len : 0.0
          end
        end
        write_region_depth_value(output, t.name, start_abs, stop_abs, nil, tid, me, mean_sum, stop_local - start_local)
        if tid != Core::CoverageResult::NoData.value
          chrom_region_stat = chrom_region_stat + Stats::DepthStat.from_array(coverage, start_local, stop_local - 1)
          representative = if @config.use_median?
                             me.to_i
                           else
                             round_nonnegative_ratio(mean_sum, stop_local - start_local)
                           end
          idx = [representative, region_dist.size - 1].min
          region_dist[idx] += 1
        end

        if @config.has_thresholds?
          thresholds = @config.threshold_values
          counts = count_threshold_bases_offset(coverage, start_abs, stop_abs, thresholds, tid, offset)
          output.write_threshold_counts(t.name, start_abs, stop_abs, nil, counts)
        end
        start_local = stop_local
      end
      chrom_region_stat
    end

    private def process_bed_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                    bed_map : Hash(String, Array(Core::Region))?,
                                    cs : Stats::IntHistogram, output : FileIO::OutputManager,
                                    region_dist : Array(Int64), offset : Int32, effective_len : Int32) : Stats::DepthStat
      chrom_region_stat = Stats::DepthStat.new
      regs = bed_map.try(&.[t.name]?) || [] of Core::Region
      region_start = offset
      region_stop = offset + effective_len
      regs.each do |region|
        s_abs = Math.max(region.start, region_start)
        e_abs = Math.min(region.stop, region_stop)
        next if e_abs <= s_abs
        s_local = s_abs - offset
        e_local = e_abs - offset

        me = 0.0
        mean_sum = 0_u64
        mean_len = e_local - s_local
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median?
            cs.clear
            (s_local...Math.min(e_local, coverage.size)).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            (s_local...Math.min(e_local, coverage.size)).each do |i|
              depth = coverage[i]
              mean_sum += depth.to_u64 if depth > 0
            end
            me = mean_len > 0 ? mean_sum.to_f / mean_len : 0.0
          end
        end
        write_region_depth_value(output, t.name, s_abs, e_abs, region.name, tid, me, mean_sum, mean_len)
        if tid != Core::CoverageResult::NoData.value && @config.window_size == 0
          chrom_region_stat = chrom_region_stat + Stats::DepthStat.from_array(coverage, s_local, e_local - 1)
          self.class.bump_distribution!(region_dist, coverage, s_local, e_local)
        end

        if @config.has_thresholds?
          thresholds = @config.threshold_values
          counts = count_threshold_bases_offset(coverage, s_abs, e_abs, thresholds, tid, offset)
          output.write_threshold_counts(t.name, s_abs, e_abs, region.name, counts)
        end
      end
      chrom_region_stat
    end

    private def round_nonnegative_ratio(sum : UInt64, length : Int32) : Int32
      return 0 if length <= 0

      numerator = sum.to_u128 * 2_u128 + length.to_u128
      rounded = numerator // (length.to_u128 * 2_u128)
      rounded > Int32::MAX ? Int32::MAX : rounded.to_i32
    end

    private def write_quantized_intervals(t : Core::Target, coverage : Core::Coverage, tid : Int32, output : FileIO::OutputManager, offset : Int32, target_size : Int32)
      quants = @config.quantize_args
      return if quants.empty?

      if tid == Core::CoverageResult::NoData.value
        # Handle case with no data - write entire chromosome as first quantize bin if it includes 0
        if quants[0] == 0
          lookup = Stats::Quantize.make_lookup(quants)
          # No data in this (sub)region
          unless lookup.empty?
            # write exactly over effective length [0, target_size-1]
            output.write_quantized_interval(t.name, offset, offset + (target_size - 1), lookup[0])
          end
        end
      else
        # Generate quantized segments using the quantize module
        Stats::Quantize.gen_quantized(quants, coverage, target_size) do |start, stop, label|
          output.write_quantized_interval(t.name, start + offset, stop + offset, label)
        end
      end
    end

    private def count_threshold_bases_offset(coverage : Core::Coverage, abs_start : Int32, abs_stop : Int32,
                                             thresholds : Array(Int32), tid : Int32, offset : Int32) : Array(Int32)
      counts = Array(Int32).new(thresholds.size, 0)
      return counts if tid == Core::CoverageResult::NoData.value
      s = (abs_start - offset).clamp(0, coverage.size)
      e = (abs_stop - offset).clamp(0, coverage.size)
      (s...e).each do |i|
        depth = coverage[i]
        thresholds.each_with_index do |threshold, idx|
          counts[idx] += 1 if depth >= threshold
        end
      end
      counts
    end
  end
end
