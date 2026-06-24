require "../stats/depth_stat"
require "../config"
require "hts"

module Depth::FileIO
  class OutputManager
    getter f_summary : File?
    getter f_global : File?
    getter f_region : File?
    getter f_perbase : (File | HTS::Bgzf)?
    getter f_regions : (File | HTS::Bgzf)?
    getter f_quantized : (File | HTS::Bgzf)?
    getter f_thresholds : (File | HTS::Bgzf)?
    # Internal line buffers for BGZF outputs
    @perbase_buf : BgzfLineBuffer?
    @regions_buf : BgzfLineBuffer?
    @quantized_buf : BgzfLineBuffer?
    @thresholds_buf : BgzfLineBuffer?
    @header_written = false
    @prefix : String
    @config : Depth::Config
    @paths_to_index : Array(String)
    @precision : Int32
    @zero_float : String
    @mean_scale : UInt128?
    @buffer_size : Int32

    # Simple aggregation buffer to reduce per-line bgzf_write overhead
    class BgzfLineBuffer
      def initialize(@bgzf : HTS::Bgzf, capacity : Int32)
        @threshold = capacity
        @buf = IO::Memory.new(capacity)
      end

      def add_line(line : String)
        @buf << line
        @buf << '\n'
        flush if @buf.size >= @threshold
      end

      def add_line(&)
        yield @buf
        @buf << '\n'
        flush if @buf.size >= @threshold
      end

      def flush
        return if @buf.size == 0
        slice = @buf.to_slice
        @bgzf.write(slice)
        @buf.clear
      end

      def close
        flush
        @bgzf.close
      end
    end

    def initialize(config : Config)
      @config = config
      @prefix = config.prefix
      label = config.mos_style? ? "mosdepth" : "mopdepth"
      @paths_to_index = [] of String
      @precision = resolve_precision
      @zero_float = zero_float_string(@precision)
      @mean_scale = decimal_scale(@precision)
      @buffer_size = (ENV["MOPDEPTH_BGZF_BUFFER"]? || "2097152").to_i

      @f_summary = File.open(path_for(label, "summary.txt"), "w")
      @f_global = File.open(path_for(label, "global.dist.txt"), "w")
      @f_region = config.has_regions? ? File.open(path_for(label, "region.dist.txt"), "w") : nil

      @f_perbase = config.no_per_base? ? nil : open_indexed_bgzf("per-base.bed.gz")
      @f_regions = config.has_regions? ? open_indexed_bgzf("regions.bed.gz") : nil
      @f_quantized = config.has_quantize? ? open_indexed_bgzf("quantized.bed.gz") : nil
      @f_thresholds = config.has_thresholds? ? open_indexed_bgzf("thresholds.bed.gz") : nil

      # Attach buffers
      @perbase_buf = wrap(@f_perbase)
      @regions_buf = wrap(@f_regions)
      @quantized_buf = wrap(@f_quantized)
      @thresholds_buf = wrap(@f_thresholds)
    end

    def write_summary_line(region : String, stat : Depth::Stats::DepthStat)
      return unless summary = @f_summary

      unless @header_written
        summary << "chrom\tlength\tbases\tmean\tmin\tmax\n"
        @header_written = true
      end

      mean = stat.n_bases > 0 ? stat.sum_depth.to_f / stat.n_bases : 0.0
      minv = stat.min_depth == Int32::MAX ? 0 : stat.min_depth
      # mosdepth uses cumulative depth in the 'bases' column
      summary << region << '\t' << stat.n_bases << '\t' << stat.sum_depth << '\t'
      write_float(summary, mean)
      summary << '\t' << minv << '\t' << stat.max_depth << '\n'
    end

    # Optionally call at end to add a total line like mosdepth
    def write_summary_total(total : Depth::Stats::DepthStat)
      return unless summary = @f_summary
      mean = total.n_bases > 0 ? total.sum_depth.to_f / total.n_bases : 0.0
      minv = total.min_depth == Int32::MAX ? 0 : total.min_depth
      summary << "total\t" << total.n_bases << '\t' << total.sum_depth << '\t'
      write_float(summary, mean)
      summary << '\t' << minv << '\t' << total.max_depth << '\n'
    end

    def write_per_base_interval(chrom : String, start : Int32, stop : Int32, depth : Int32)
      return unless perbase = @f_perbase
      if buf = @perbase_buf
        buf.add_line do |io|
          io << chrom << '\t' << start << '\t' << stop << '\t' << depth
        end
      else
        write_raw_line(perbase) do |io|
          io << chrom << '\t' << start << '\t' << stop << '\t' << depth
        end
      end
    end

    def write_region_stat(chrom : String, start : Int32, stop : Int32, name : String?, value : Float64)
      write_region_value(chrom, start, stop, name) do |io|
        write_float(io, value)
      end
    end

    def write_region_mean(chrom : String, start : Int32, stop : Int32, name : String?, sum : UInt64, length : Int32)
      write_region_value(chrom, start, stop, name) do |io|
        write_mean(io, sum, length)
      end
    end

    def write_region_zero(chrom : String, start : Int32, stop : Int32, name : String?)
      write_region_value(chrom, start, stop, name) do |io|
        io << @zero_float
      end
    end

    private def write_region_value(chrom : String, start : Int32, stop : Int32, name : String?, &)
      return unless regions = @f_regions

      if buf = @regions_buf
        buf.add_line do |io|
          write_region_fields(io, chrom, start, stop, name)
          yield io
        end
      else
        write_raw_line(regions) do |io|
          write_region_fields(io, chrom, start, stop, name)
          yield io
        end
      end
    end

    def write_quantized_interval(chrom : String, start : Int32, stop : Int32, label : String)
      return unless quantized = @f_quantized
      if buf = @quantized_buf
        buf.add_line do |io|
          io << chrom << '\t' << start << '\t' << stop << '\t' << label
        end
      else
        write_raw_line(quantized) do |io|
          io << chrom << '\t' << start << '\t' << stop << '\t' << label
        end
      end
    end

    def write_thresholds_header(thresholds : Array(Int32))
      return unless thresholds_io = @f_thresholds
      if buf = @thresholds_buf
        buf.add_line do |io|
          write_threshold_header_fields(io, thresholds)
        end
      else
        write_raw_line(thresholds_io) do |io|
          write_threshold_header_fields(io, thresholds)
        end
      end
    end

    def write_threshold_counts(chrom : String, start : Int32, stop : Int32, name : String?, counts : Array(Int32))
      return unless thresholds_io = @f_thresholds
      if buf = @thresholds_buf
        buf.add_line do |io|
          write_threshold_count_fields(io, chrom, start, stop, name, counts)
        end
      else
        write_raw_line(thresholds_io) do |io|
          write_threshold_count_fields(io, chrom, start, stop, name, counts)
        end
      end
    end

    def close_all
      [@f_summary, @f_global, @f_region].each(&.try(&.close))
      # Flush and close buffered BGZF streams
      @perbase_buf.try(&.close)
      @regions_buf.try(&.close)
      @quantized_buf.try(&.close)
      @thresholds_buf.try(&.close)
      # Close any BGZF streams not wrapped in a buffer (fallback)
      [@f_perbase, @f_regions, @f_quantized, @f_thresholds].each do |io|
        io.try do |file|
          begin
            file.close
          rescue
          end
        end
      end

      # Always build CSI indices for BGZF interval outputs
      build_csi_indices
    end

    private def build_csi_indices
      # Use htslib's tbx_index_build3 with a BED preset configuration.
      # To avoid the issue where references to $tbx_conf_bed become __imp_* with MinGW static linking,
      # construct TbxConfT locally and pass it instead of using a global variable.
      # min_shift=14 requests CSI index (instead of TBI)
      bed_conf = HTS::LibHTS::TbxConfT.new
      bed_conf.preset = 0x10000    # TBX_UCSC (BED-like 0-based)
      bed_conf.sc = 1              # seq col (1-based)
      bed_conf.bc = 2              # begin col (1-based)
      bed_conf.ec = 3              # end col (1-based)
      bed_conf.meta_char = '#'.ord # comment/meta line start
      bed_conf.line_skip = 0
      conf_ptr = pointerof(bed_conf)

      @paths_to_index.each do |gz_path|
        # Build explicit .csi next to gz
        csi = "#{gz_path}.csi"
        ret = HTS::LibHTS.tbx_index_build3(gz_path, csi, 14, @config.threads, conf_ptr)
        if ret != 0
          STDERR.puts "[mopdepth] warning: failed to build CSI for #{gz_path} (tbx_index_build3 ret=#{ret})"
        end
      end
    end

    # Helper methods
    private def path_for(label : String, suffix : String) : String
      "#{@prefix}.#{label}.#{suffix}"
    end

    private def open_indexed_bgzf(basename : String)
      path = "#{@prefix}.#{basename}"
      @paths_to_index << path
      HTS::Bgzf.open(path, "w1z", threads: @config.threads)
    end

    private def wrap(io : (File | HTS::Bgzf)?)
      return nil unless io.is_a?(HTS::Bgzf)
      BgzfLineBuffer.new(io.as(HTS::Bgzf), @buffer_size)
    end

    private def write_float(io, value : Float64)
      if value == 0.0
        io << @zero_float
        return
      end

      value.format(io, decimal_places: @precision)
    end

    private def write_mean(io, sum : UInt64, length : Int32)
      if sum == 0 || length <= 0
        io << @zero_float
        return
      end

      unless scale = @mean_scale
        write_float(io, sum.to_f / length)
        return
      end

      len = length.to_u128
      scaled = (sum.to_u128 * scale + (len // 2)) // len
      whole = scaled // scale
      fraction = scaled % scale

      io << whole
      return if @precision <= 0

      io << '.'
      write_padded_fraction(io, fraction, @precision)
    end

    private def write_raw_line(io : File, &)
      yield io
      io << '\n'
    end

    private def write_raw_line(io : HTS::Bgzf, &)
      line = String.build do |s|
        yield s
      end
      io.puts(line)
    end

    private def write_region_fields(io, chrom : String, start : Int32, stop : Int32, name : String?)
      io << chrom << '\t' << start << '\t' << stop << '\t'
      io << name << '\t' if name
    end

    private def write_threshold_header_fields(io, thresholds : Array(Int32))
      io << "#chrom\tstart\tend\tregion"
      thresholds.each { |threshold| io << '\t' << threshold << 'X' }
    end

    private def write_threshold_count_fields(io, chrom : String, start : Int32, stop : Int32,
                                             name : String?, counts : Array(Int32))
      io << chrom << '\t' << start << '\t' << stop << '\t' << (name || "unknown")
      counts.each { |count| io << '\t' << count }
    end

    private def resolve_precision : Int32
      (ENV["MOPDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
    end

    private def zero_float_string(precision : Int32) : String
      return "0" if precision <= 0

      String.build do |io|
        io << "0."
        precision.times { io << '0' }
      end
    end

    private def decimal_scale(precision : Int32) : UInt128?
      return 1_u128 if precision <= 0
      return nil if precision > 18

      scale = 1_u128
      precision.times { scale *= 10_u128 }
      scale
    end

    private def write_padded_fraction(io, fraction : UInt128, precision : Int32)
      divisor = 1_u128
      (precision - 1).times { divisor *= 10_u128 }

      while divisor > 0
        digit = fraction // divisor
        io << digit.to_i
        fraction %= divisor
        divisor //= 10_u128
      end
    end
  end
end
