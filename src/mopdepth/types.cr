module Mopdepth
  enum Mode
    Standard
    Fast
    Fragment
  end

  # Filtering and counting behavior for coverage calculations.
  #
  # Fragment lengths are absolute template lengths. `nil` means unbounded.
  record Options,
    mode : Mode = Mode::Standard,
    minimum_mapping_quality : Int32 = 0,
    minimum_fragment_length : Int32? = nil,
    maximum_fragment_length : Int32? = nil,
    exclude_flags : UInt16 = 1796_u16,
    include_flags : UInt16 = 0_u16,
    read_groups : Array(String) = [] of String do
    def validate! : Nil
      unless (0..255).includes?(minimum_mapping_quality)
        raise ArgumentError.new("minimum_mapping_quality must be between 0 and 255")
      end
      if minimum = minimum_fragment_length
        raise ArgumentError.new("minimum_fragment_length must not be negative") if minimum < 0
      end
      if maximum = maximum_fragment_length
        raise ArgumentError.new("maximum_fragment_length must not be negative") if maximum < 0
      end
      if (minimum = minimum_fragment_length) && (maximum = maximum_fragment_length) && minimum > maximum
        raise ArgumentError.new("minimum_fragment_length must not exceed maximum_fragment_length")
      end
    end
  end

  record Target, name : String, length : Int32, tid : Int32 do
    def region : Region
      Region.new(name, 0, length)
    end
  end

  # A zero-based, half-open genomic interval.
  record Region, chromosome : String, start : Int32, stop : Int32 do
    def initialize(@chromosome : String, @start : Int32, @stop : Int32)
      raise ArgumentError.new("chromosome must not be empty") if chromosome.empty?
      raise ArgumentError.new("region start must not be negative") if start < 0
      raise ArgumentError.new("region stop must not precede start") if stop < start
    end

    def length : Int32
      stop - start
    end

    def empty? : Bool
      start == stop
    end
  end

  enum CoverageStatus
    Observed
    NoData
  end

  # Dense per-base coverage for `region`. `values` is caller-owned when the
  # result came from `Reader#coverage_into`.
  record Coverage,
    region : Region,
    values : Array(Int32),
    status : CoverageStatus,
    records_examined : Int64 do
    def observed? : Bool
      status.observed?
    end
  end

  enum Stage
    Preparing
    Reading
    Finalizing
  end

  record Progress,
    stage : Stage,
    region : Region? = nil,
    target_index : Int32 = 0,
    target_count : Int32 = 0,
    records_examined : Int64 = 0_i64
end
