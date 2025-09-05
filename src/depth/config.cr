require "./stats/quantize"
require "./stats/threshold"
require "./core/options"

module Depth
  class Config
    property prefix : String = ""
    property path : String = ""
    property threads : Int32 = 0
    property chrom : String = ""
    property by : String = "" # numeric window or BED path
    property? no_per_base : Bool = false
    property mapq : Int32 = 0
    property min_frag_len : Int32 = -1
    property max_frag_len : Int32 = -1
    property? fast_mode : Bool = false
    property? fragment_mode : Bool = false
    property? use_median : Bool = false
    property thresholds : Array(Int32) = [] of Int32
    property thresholds_str : String = ""
    property quantize : String = ""
    property exclude_flag : UInt16 = 1796_u16
    property include_flag : UInt16 = 0_u16
    property read_groups_str : String = ""
    # When true, use mosdepth-compatible names (mosdepth.*); otherwise mopdepth.*
    property? mos_style : Bool = false

    def validate!
      raise ArgumentError.new("BAM/CRAM path is required") if path.empty?
      raise ArgumentError.new("Output prefix is required") if prefix.empty?
      raise ArgumentError.new("Invalid MAPQ threshold") if mapq < 0
      raise ArgumentError.new("Invalid thread count") if threads < 0

      if min_frag_len >= 0 && max_frag_len >= 0 && min_frag_len > max_frag_len
        raise ArgumentError.new("min_frag_len cannot be greater than max_frag_len")
      end

      if fast_mode? && fragment_mode?
        raise ArgumentError.new("--fast-mode and --fragment-mode cannot be used together")
      end

      if has_thresholds? && !has_regions?
        raise ArgumentError.new("--thresholds can only be used when --by is specified")
      end
    end

    def to_options : Core::Options
      Core::Options.new(
        mapq: mapq,
        min_frag_len: min_frag_len,
        max_frag_len: (max_frag_len < 0 ? Int32::MAX : max_frag_len),
        exclude_flag: exclude_flag,
        include_flag: include_flag,
        fast_mode: fast_mode?,
        fragment_mode: fragment_mode?,
        read_groups: parse_read_groups,
      )
    end

    def has_regions? : Bool
      !by.empty?
    end

    def window_size : Int32
      return 0 unless has_regions?
      return 0 unless by.each_char.all?(&.ascii_number?)
      by.to_i
    end

    def bed_path : String?
      return nil unless has_regions?
      return nil if by.each_char.all?(&.ascii_number?)
      by
    end

    def has_quantize? : Bool
      !quantize.empty? && quantize != "nil"
    end

    def quantize_args : Array(Int32)
      return [] of Int32 unless has_quantize?
      Stats::Quantize.get_quantize_args(quantize)
    end

    def has_thresholds? : Bool
      !thresholds_str.empty? && thresholds_str != "nil"
    end

    def threshold_values : Array(Int32)
      return [] of Int32 unless has_thresholds?
      Stats.threshold_args(thresholds_str)
    end

    private def parse_read_groups : Array(String)
      return [] of String if read_groups_str.empty? || read_groups_str == "nil"
      read_groups_str.split(',').map(&.strip)
    end
  end
end
