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
    @header_written = false
    @prefix : String
    @config : Depth::Config
    @paths_to_index : Array(String)

    def initialize(config : Config)
      @config = config
      @prefix = config.prefix
      label = config.mos_style? ? "mosdepth" : "mopdepth"

      @f_summary = File.open("#{@prefix}.#{label}.summary.txt", "w")
      @f_global = File.open("#{@prefix}.#{label}.global.dist.txt", "w")
      @f_region = config.has_regions? ? File.open("#{@prefix}.#{label}.region.dist.txt", "w") : nil
      # Use BGZF for BED-like interval outputs (mosdepth-compatible behavior)
      @paths_to_index = [] of String
      @f_perbase = if config.no_per_base?
                     nil
                   else
                     path = "#{@prefix}.per-base.bed.gz"
                     @paths_to_index << path
                     HTS::Bgzf.open(path, "wz")
                   end

      @f_regions = if config.has_regions?
                     path = "#{@prefix}.regions.bed.gz"
                     @paths_to_index << path
                     HTS::Bgzf.open(path, "wz")
                   else
                     nil
                   end

      @f_quantized = if config.has_quantize?
                       path = "#{@prefix}.quantized.bed.gz"
                       @paths_to_index << path
                       HTS::Bgzf.open(path, "wz")
                     else
                       nil
                     end

      @f_thresholds = if config.has_thresholds?
                        path = "#{@prefix}.thresholds.bed.gz"
                        @paths_to_index << path
                        HTS::Bgzf.open(path, "wz")
                      else
                        nil
                      end
    end

    def write_summary_line(region : String, stat : Depth::Stats::DepthStat)
      return unless @f_summary

      unless @header_written
        @f_summary.not_nil! << ["chrom", "length", "bases", "mean", "min", "max"].join("\t") << '\n'
        @header_written = true
      end

      mean = stat.n_bases > 0 ? stat.sum_depth.to_f / stat.n_bases : 0.0
      pr = (ENV["MOPDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      mean_str = sprintf("%.#{pr}f", mean)
      minv = stat.min_depth == Int32::MAX ? 0 : stat.min_depth
      # mosdepth uses cumulative depth in the 'bases' column
      @f_summary.not_nil! << [region, stat.n_bases, stat.sum_depth, mean_str, minv, stat.max_depth].join("\t") << '\n'
    end

    # Optionally call at end to add a total line like mosdepth
    def write_summary_total(total : Depth::Stats::DepthStat)
      return unless @f_summary
      mean = total.n_bases > 0 ? total.sum_depth.to_f / total.n_bases : 0.0
      pr = (ENV["MOPDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      mean_str = sprintf("%.#{pr}f", mean)
      minv = total.min_depth == Int32::MAX ? 0 : total.min_depth
      @f_summary.not_nil! << ["total", total.n_bases, total.sum_depth, mean_str, minv, total.max_depth].join("\t") << "\n"
    end

    def write_per_base_interval(chrom : String, start : Int32, stop : Int32, depth : Int32)
      return unless @f_perbase
      @f_perbase.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{depth}")
    end

    def write_region_stat(chrom : String, start : Int32, stop : Int32, name : String?, value : Float64)
      return unless @f_regions
      pr = (ENV["MOPDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      val_str = sprintf("%.#{pr}f", value)
      if name
        @f_regions.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{name}\t#{val_str}")
      else
        @f_regions.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{val_str}")
      end
    end

    def write_quantized_interval(chrom : String, start : Int32, stop : Int32, label : String)
      return unless @f_quantized
      @f_quantized.not_nil!.puts("#{chrom}\t#{start}\t#{stop}\t#{label}")
    end

    def write_thresholds_header(thresholds : Array(Int32))
      return unless @f_thresholds
      line = String.build do |io|
        io << "#chrom\tstart\tend\tregion"
        thresholds.each { |t| io << "\t" << t << "X" }
      end
      @f_thresholds.not_nil!.puts(line)
    end

    def write_threshold_counts(chrom : String, start : Int32, stop : Int32, name : String?, counts : Array(Int32))
      return unless @f_thresholds
      line = String.build do |io|
        io << chrom << '\t' << start << '\t' << stop << '\t' << (name || "unknown")
        counts.each { |c| io << '\t' << c }
      end
      @f_thresholds.not_nil!.puts(line)
    end

    def close_all
      [@f_summary, @f_global, @f_region].each(&.try(&.close))
      # BGZF streams (respond to close as well)
      [@f_perbase, @f_regions, @f_quantized, @f_thresholds].each do |io|
        io.try(&.close)
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
  end
end
