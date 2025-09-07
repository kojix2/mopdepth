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

      def flush
        return if @buf.size == 0
        slice = @buf.to_slice
        @bgzf.write(slice)
        @buf = IO::Memory.new(@threshold)
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
      return unless @f_summary

      unless @header_written
        @f_summary.not_nil! << ["chrom", "length", "bases", "mean", "min", "max"].join("\t") << '\n'
        @header_written = true
      end

      mean = stat.n_bases > 0 ? stat.sum_depth.to_f / stat.n_bases : 0.0
      mean_str = sprintf("%.#{@precision}f", mean)
      minv = stat.min_depth == Int32::MAX ? 0 : stat.min_depth
      # mosdepth uses cumulative depth in the 'bases' column
      @f_summary.not_nil! << [region, stat.n_bases, stat.sum_depth, mean_str, minv, stat.max_depth].join("\t") << '\n'
    end

    # Optionally call at end to add a total line like mosdepth
    def write_summary_total(total : Depth::Stats::DepthStat)
      return unless @f_summary
      mean = total.n_bases > 0 ? total.sum_depth.to_f / total.n_bases : 0.0
      mean_str = sprintf("%.#{@precision}f", mean)
      minv = total.min_depth == Int32::MAX ? 0 : total.min_depth
      @f_summary.not_nil! << ["total", total.n_bases, total.sum_depth, mean_str, minv, total.max_depth].join("\t") << "\n"
    end

    def write_per_base_interval(chrom : String, start : Int32, stop : Int32, depth : Int32)
      return unless @f_perbase
      if buf = @perbase_buf
        buf.add_line("#{chrom}\t#{start}\t#{stop}\t#{depth}")
      else
        @f_perbase.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{depth}")
      end
    end

    def write_region_stat(chrom : String, start : Int32, stop : Int32, name : String?, value : Float64)
      return unless @f_regions
      val_str = sprintf("%.#{@precision}f", value)
      if name
        if buf = @regions_buf
          buf.add_line("#{chrom}\t#{start}\t#{stop}\t#{name}\t#{val_str}")
        else
          @f_regions.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{name}\t#{val_str}")
        end
      else
        if buf = @regions_buf
          buf.add_line("#{chrom}\t#{start}\t#{stop}\t#{val_str}")
        else
          @f_regions.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{val_str}")
        end
      end
    end

    def write_quantized_interval(chrom : String, start : Int32, stop : Int32, label : String)
      return unless @f_quantized
      if buf = @quantized_buf
        buf.add_line("#{chrom}\t#{start}\t#{stop}\t#{label}")
      else
        @f_quantized.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{label}")
      end
    end

    def write_thresholds_header(thresholds : Array(Int32))
      return unless @f_thresholds
      line = String.build do |io|
        io << "#chrom\tstart\tend\tregion"
        thresholds.each { |t| io << "\t" << t << "X" }
      end
      if buf = @thresholds_buf
        buf.add_line(line)
      else
        @f_thresholds.not_nil!.puts(line)
      end
    end

    def write_threshold_counts(chrom : String, start : Int32, stop : Int32, name : String?, counts : Array(Int32))
      return unless @f_thresholds
      line = String.build do |io|
        io << chrom << '\t' << start << '\t' << stop << '\t' << (name || "unknown")
        counts.each { |c| io << '\t' << c }
      end
      if buf = @thresholds_buf
        buf.add_line(line)
      else
        @f_thresholds.not_nil!.puts(line)
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
        io.try do |f|
          begin
            f.close
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

      @paths_to_index.each do |gz|
        # Build explicit .csi next to gz
        csi = "#{gz}.csi"
        ret = HTS::LibHTS.tbx_index_build3(gz, csi, 14, @config.threads, conf_ptr)
        if ret != 0
          STDERR.puts "[mopdepth] warning: failed to build CSI for #{gz} (tbx_index_build3 ret=#{ret})"
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
      HTS::Bgzf.open(path, "w1z")
    end

    private def wrap(io : (File | HTS::Bgzf)?)
      return nil unless io.is_a?(HTS::Bgzf)
      BgzfLineBuffer.new(io.as(HTS::Bgzf), @buffer_size)
    end

    private def resolve_precision : Int32
      (ENV["MOPDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
    end
  end
end
