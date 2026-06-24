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

      # Open BAM/CRAM
      bam = HTS::Bam.open(@config.path, threads: @config.threads)
      begin
        bam.load_index
      rescue ex
        raise "Failed to load index for #{@config.path}: #{ex.message}"
      end

      region = FileIO.parse_region_str(@config.chrom)
      opts = @config.to_options

      # Create output manager
      output = FileIO::OutputManager.new(@config)

      # Write threshold header if needed
      if @config.has_thresholds?
        thresholds = @config.threshold_values
        output.write_thresholds_header(thresholds)
      end

      begin
        # Get reference sequences from header using hts.cr API
        target_names = bam.header.target_names
        target_lengths = bam.header.target_len
        targets = target_names.map_with_index do |name, i|
          Core::Target.new(name, target_lengths[i].to_i32, i)
        end

        sub_targets = if selected_region = region
                        selected = targets.select { |target| target.name == selected_region.chrom }
                        # If a chromosome was specified but not found, fail (tests expect non-zero status)
                        if selected.empty?
                          raise ConfigError.new("Chromosome not found: #{selected_region.chrom}")
                        end
                        selected
                      else
                        targets
                      end

        # Initialize statistics
        global_dist = Array(Int64).new(512, 0_i64)
        total_global_dist = Array(Int64).new(512, 0_i64)
        region_dist = Array(Int64).new(512, 0_i64)
        total_region_dist = Array(Int64).new(512, 0_i64)
        global_stat = Stats::DepthStat.new
        global_region_stat = Stats::DepthStat.new
        cs = Stats::IntHistogram.new(@config.use_median? ? 65_536 : 0)

        # Handle window/BED regions
        bed_map : Hash(String, Array(Core::Region))? = nil
        window = @config.window_size
        if bed_path = @config.bed_path
          bed_map = FileIO.read_bed(bed_path)
        end

        # Create coverage calculator
        calculator = Core::CoverageCalculator.new(bam, opts)

        # Reusable coverage buffer (grow-only). We'll only use [0, effective_len]
        coverage = Core::Coverage.new(0)
        coverage_dirty = false

        # Process each target
        sub_targets.each do |target|
          # Skip if no regions for this chrom and we won't write per-base
          if @config.no_per_base? && (regions_by_chrom = bed_map) && !regions_by_chrom.has_key?(target.name)
            next
          end

          # Determine query region and sizing
          query_region = if selected_region = region
                           target.name == selected_region.chrom ? selected_region : Core::Region.new(target.name, 0, 0)
                         else
                           Core::Region.new(target.name, 0, 0)
                         end

          # If a region is provided, shrink coverage to [start, stop); otherwise use full chromosome
          offset = 0
          effective_len = target.length
          if region && target.name == query_region.chrom && (query_region.start > 0 || query_region.stop > 0)
            offset = query_region.start
            stop = (query_region.stop > 0 ? query_region.stop : target.length)
            effective_len = (stop - offset)
          end

          target_size = effective_len + 1
          if coverage.size < target_size
            # grow once; keep capacity for reuse
            coverage.concat(Array(Int32).new(target_size - coverage.size, 0))
          end
          if coverage_dirty
            i_full = 0
            while i_full < target_size
              coverage[i_full] = 0
              i_full += 1
            end
            coverage_dirty = false
          end

          tid = calculator.calculate(coverage, query_region, offset)
          next if tid == Core::CoverageResult::ChromNotFound.value

          # Build final coverage from diff-array for [0, effective_len]
          if tid != Core::CoverageResult::NoData.value
            i = 0
            sum = 0
            while i < target_size
              sum += coverage[i]
              coverage[i] = sum
              i += 1
            end
            coverage_dirty = true
          end

          # Write per-base intervals
          if output.f_perbase
            if tid == Core::CoverageResult::NoData.value
              write_len = (region && target.name == query_region.chrom) ? effective_len : target.length
              output.write_per_base_interval(target.name, offset, offset + write_len, 0)
            else
              # restrict per-base segments to [0, effective_len]
              self.class.each_constant_segment(coverage, target_size - 1) do |(s, e, v)|
                output.write_per_base_interval(target.name, s + offset, e + offset, v)
              end
            end
          end

          # Write quantized intervals
          if output.f_quantized && @config.has_quantize?
            write_quantized_intervals(target, coverage, tid, output, offset, target_size)
          end

          # Process regions (window or BED)
          chrom_region_stat = Stats::DepthStat.new
          if output.f_regions
            chrom_region_stat = write_region_stats_with_offset(target, coverage, tid, window, bed_map, cs, output, region_dist, offset, effective_len)
          end

          # Process per-chromosome distributions and stats
          if tid != Core::CoverageResult::NoData.value
            self.class.bump_distribution!(global_dist, coverage, 0, target_size - 1)
            chrom_stat = Stats::DepthStat.from_array(coverage, 0, target_size - 2)
            global_stat = global_stat + chrom_stat
            output.write_summary_line(target.name, chrom_stat)
            if output.f_regions
              global_region_stat = global_region_stat + chrom_region_stat
              output.write_summary_line("#{target.name}_region", chrom_region_stat)
            end
          end

          # Write distributions
          if f_global = output.f_global
            self.class.write_distribution(f_global.as(::IO), target.name, global_dist)
            # accumulate into genome-wide total
            self.class.sum_into!(total_global_dist, global_dist)
          end
          if f_region = output.f_region
            self.class.write_distribution(f_region.as(::IO), target.name, region_dist)
            # accumulate into genome-wide total for regions
            self.class.sum_into!(total_region_dist, region_dist)
            region_dist.fill(0_i64)
          end
          global_dist.fill(0_i64)
        end
        # Append mosdepth-like total lines
        output.write_summary_total(global_stat)
        output.write_summary_line("total_region", global_region_stat) if output.f_regions
        if f_global = output.f_global
          self.class.write_distribution(f_global.as(::IO), "total", total_global_dist)
        end
        if f_region = output.f_region
          self.class.write_distribution(f_region.as(::IO), "total", total_region_dist)
        end
      ensure
        output.close_all
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
