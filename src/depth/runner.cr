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
      @config.validate!
      region = FileIO.parse_region_str(@config.chrom)
      opts = @config.to_options
      bam = open_bam

      output = FileIO::OutputManager.new(@config)
      write_threshold_header(output)

      begin
        targets = selected_targets(load_targets(bam), region)
        bed_map = load_bed_map
        state = ProcessingState.new(@config.use_median?)
        process_targets(bam, opts, targets, region, bed_map, output, state)
        write_total_outputs(output, state)
      ensure
        output.close_all
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

    private def open_bam : HTS::Bam
      bam = HTS::Bam.open(@config.path, threads: @config.threads)
      begin
        bam.load_index
      rescue ex
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
        Core::Target.new(name, target_lengths[i].to_i32, i)
      end
    end

    private def selected_targets(targets : Array(Core::Target), region : Core::Region?) : Array(Core::Target)
      return targets unless selected_region = region

      selected = targets.select { |target| target.name == selected_region.chrom }
      raise ConfigError.new("Chromosome not found: #{selected_region.chrom}") if selected.empty?
      selected
    end

    private def load_bed_map : Hash(String, Array(Core::Region))?
      return unless bed_path = @config.bed_path

      FileIO.read_bed(bed_path)
    end

    private def process_targets(bam : HTS::Bam, opts : Core::Options, targets : Array(Core::Target),
                                region : Core::Region?, bed_map : Hash(String, Array(Core::Region))?,
                                output : FileIO::OutputManager, state : ProcessingState)
      calculator = Core::CoverageCalculator.new(bam, opts)
      coverage = Core::Coverage.new(0)
      coverage_dirty = false
      window = @config.window_size

      targets.each do |target|
        next if skip_target?(target, bed_map)

        slice = target_slice(target, region)
        prepare_coverage!(coverage, slice[:target_size], coverage_dirty)
        tid = calculator.calculate(coverage, slice[:query_region], slice[:offset])
        next if tid == Core::CoverageResult::ChromNotFound.value

        coverage_dirty = finalize_coverage!(coverage, tid, slice[:target_size])
        process_target_outputs(target, coverage, tid, slice, window, bed_map, output, state)
      end
    end

    private def skip_target?(target : Core::Target, bed_map : Hash(String, Array(Core::Region))?) : Bool
      return false unless @config.no_per_base? && bed_map

      !bed_map.has_key?(target.name)
    end

    private def target_slice(target : Core::Target, region : Core::Region?)
      query_region = if selected_region = region
                       target.name == selected_region.chrom ? selected_region : Core::Region.new(target.name, 0, 0)
                     else
                       Core::Region.new(target.name, 0, 0)
                     end

      offset = 0
      effective_len = target.length
      if region && target.name == query_region.chrom && (query_region.start > 0 || query_region.stop > 0)
        offset = query_region.start
        stop = (query_region.stop > 0 ? query_region.stop : target.length)
        effective_len = stop - offset
      end

      {query_region: query_region, offset: offset, effective_len: effective_len, target_size: effective_len + 1}
    end

    private def prepare_coverage!(coverage : Core::Coverage, target_size : Int32, coverage_dirty : Bool)
      if coverage.size < target_size
        coverage.concat(Array(Int32).new(target_size - coverage.size, 0))
      end
      clear_coverage!(coverage, target_size) if coverage_dirty
    end

    private def clear_coverage!(coverage : Core::Coverage, target_size : Int32)
      i_full = 0
      while i_full < target_size
        coverage[i_full] = 0
        i_full += 1
      end
    end

    private def finalize_coverage!(coverage : Core::Coverage, tid : Int32, target_size : Int32) : Bool
      return false if tid == Core::CoverageResult::NoData.value

      i = 0
      sum = 0
      while i < target_size
        sum += coverage[i]
        coverage[i] = sum
        i += 1
      end
      true
    end

    private def process_target_outputs(target : Core::Target, coverage : Core::Coverage, tid : Int32, slice,
                                       window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                       output : FileIO::OutputManager, state : ProcessingState)
      write_per_base_intervals(target, coverage, tid, slice, output)
      write_quantized_intervals(target, coverage, tid, output, slice[:offset], slice[:target_size]) if output.f_quantized && @config.has_quantize?

      chrom_region_stat = write_regions_if_needed(target, coverage, tid, slice, window, bed_map, output, state)
      write_target_stats(target, coverage, tid, slice[:target_size], chrom_region_stat, output, state)
      write_target_distributions(target, output, state)
    end

    private def write_per_base_intervals(target : Core::Target, coverage : Core::Coverage, tid : Int32,
                                         slice, output : FileIO::OutputManager)
      return unless output.f_perbase

      if tid == Core::CoverageResult::NoData.value
        write_len = @config.chrom.empty? ? target.length : slice[:effective_len]
        output.write_per_base_interval(target.name, slice[:offset], slice[:offset] + write_len, 0)
      else
        self.class.each_constant_segment(coverage, slice[:target_size] - 1) do |(s, e, v)|
          output.write_per_base_interval(target.name, s + slice[:offset], e + slice[:offset], v)
        end
      end
    end

    private def write_regions_if_needed(target : Core::Target, coverage : Core::Coverage, tid : Int32, slice,
                                        window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                        output : FileIO::OutputManager, state : ProcessingState) : Stats::DepthStat
      return Stats::DepthStat.new unless output.f_regions

      write_region_stats_with_offset(target, coverage, tid, window, bed_map, state.cs, output,
        state.region_dist, slice[:offset], slice[:effective_len])
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

    private def bump_threshold_counts!(counts : Array(Int32), thresholds : Array(Int32), depth : Int32, len : Int32)
      return if len <= 0 || thresholds.empty?
      thresholds.each_with_index do |threshold, idx|
        break if depth < threshold
        counts[idx] += len
      end
    end

    private def bump_depth_count!(dist : Array(Int64), depth : Int32, len : Int32)
      return if depth < 0 || len <= 0
      v = depth > Core::MAX_COVERAGE ? Core::MAX_COVERAGE - 10 : depth
      if v >= dist.size
        old = dist.size
        new_size = v + 10
        dist.concat(Array.new(new_size - old, 0_i64))
      end
      dist[v] += len
    end

    private def write_region_stats(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                   window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                   cs : Stats::IntHistogram, output : FileIO::OutputManager,
                                   region_dist : Array(Int64))
      if window > 0
        process_window_regions(t, coverage, tid, window, cs, output, region_dist)
      else
        process_bed_regions(t, coverage, tid, bed_map, cs, output, region_dist)
      end
    end

    # Offset-aware variant for region-shrunk arrays
    private def write_region_stats_with_offset(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                               window : Int32, bed_map : Hash(String, Array(Core::Region))?,
                                               cs : Stats::IntHistogram, output : FileIO::OutputManager,
                                               region_dist : Array(Int64), offset : Int32, effective_len : Int32) : Stats::DepthStat
      if window > 0
        process_window_regions_with_offset(t, coverage, tid, window, cs, output, region_dist, offset, effective_len)
      else
        process_bed_regions_with_offset(t, coverage, tid, bed_map, cs, output, region_dist, offset, effective_len)
      end
    end

    private def process_window_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                       window : Int32, cs : Stats::IntHistogram,
                                       output : FileIO::OutputManager, region_dist : Array(Int64))
      start = 0
      while start < t.length
        stop = Math.min(start + window, t.length)
        me = 0.0
        mean_sum = 0_u64
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median?
            cs.clear
            (start...stop).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            len = stop - start
            (start...stop).each do |i|
              depth = coverage[i]
              mean_sum += depth.to_u64 if depth > 0
            end
            me = len > 0 ? mean_sum.to_f / len : 0.0
          end
        end
        write_region_depth_value(output, t.name, start, stop, nil, tid, me, mean_sum, stop - start)
        if tid != Core::CoverageResult::NoData.value
          idx = [me.to_i, region_dist.size - 1].min
          region_dist[idx] += 1
        end

        # Process thresholds for this window
        if @config.has_thresholds?
          thresholds = @config.threshold_values
          counts = count_threshold_bases(coverage, start, stop, thresholds, tid)
          output.write_threshold_counts(t.name, start, stop, nil, counts)
        end
        start = stop
      end
    end

    private def process_window_regions_with_offset(t : Core::Target, coverage : Core::Coverage, tid : Int32,
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
          idx = [me.round.to_i, region_dist.size - 1].min
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
                                    region_dist : Array(Int64))
      regs = bed_map.try(&.[t.name]?) || [] of Core::Region
      regs.each do |region|
        me = 0.0
        mean_sum = 0_u64
        mean_len = region.stop - region.start
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median?
            cs.clear
            (region.start...Math.min(region.stop, coverage.size)).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            (region.start...Math.min(region.stop, coverage.size)).each do |i|
              depth = coverage[i]
              mean_sum += depth.to_u64 if depth > 0
            end
            me = mean_len > 0 ? mean_sum.to_f / mean_len : 0.0
          end
        end
        write_region_depth_value(output, t.name, region.start, region.stop, region.name, tid, me, mean_sum, mean_len)
        if tid != Core::CoverageResult::NoData.value && @config.window_size == 0
          self.class.bump_distribution!(region_dist, coverage, region.start, region.stop)
        end

        # Process thresholds for this BED region
        if @config.has_thresholds?
          thresholds = @config.threshold_values
          counts = count_threshold_bases(coverage, region.start, region.stop, thresholds, tid)
          output.write_threshold_counts(t.name, region.start, region.stop, region.name, counts)
        end
      end
    end

    private def process_bed_regions_with_offset(t : Core::Target, coverage : Core::Coverage, tid : Int32,
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

    private def count_threshold_bases(coverage : Core::Coverage, start : Int32, stop : Int32,
                                      thresholds : Array(Int32), tid : Int32) : Array(Int32)
      counts = Array(Int32).new(thresholds.size, 0)

      if tid == Core::CoverageResult::NoData.value
        # No data case - all counts are 0
        return counts
      end

      # Count bases that meet each threshold
      (start...Math.min(stop, coverage.size)).each do |i|
        depth = coverage[i]
        thresholds.each_with_index do |threshold, idx|
          counts[idx] += 1 if depth >= threshold
        end
      end

      counts
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
