module Depth::Core
  # The containing BED map owns the chromosome name once per reference.
  record Interval, start : Int32, stop : Int32, name : String? = nil
end
