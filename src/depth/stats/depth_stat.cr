module Depth::Stats
  struct DepthStat
    property n_bases : Int32 = 0
    property sum_depth : Int64 = 0_i64
    property min_depth : Int32 = Int32::MAX
    property max_depth : Int32 = 0
    # number of bases with coverage > 0 (mosdepth "bases")
    property bases : Int32 = 0

    def clear
      @n_bases = 0
      @sum_depth = 0
      @min_depth = Int32::MAX
      @max_depth = 0
      @bases = 0
    end

    def self.from_slice(slice : Slice(Int32))
      s = DepthStat.new
      s.n_bases = slice.size
      slice.each do |v|
        s.sum_depth += v
        s.min_depth = v if v < s.min_depth
        s.max_depth = v if v > s.max_depth
        s.bases += 1 if v > 0
      end
      s
    end

    def self.from_array(array : Array(Int32), start_idx : Int32 = 0, end_idx : Int32 = -1)
      s = DepthStat.new
      end_idx = array.size - 1 if end_idx < 0
      s.n_bases = end_idx - start_idx + 1
      (start_idx..end_idx).each do |i|
        v = array[i]
        s.sum_depth += v
        s.min_depth = v if v < s.min_depth
        s.max_depth = v if v > s.max_depth
        s.bases += 1 if v > 0
      end
      s
    end

    def +(other : DepthStat) : DepthStat
      result = DepthStat.new
      result.n_bases = self.n_bases + other.n_bases
      result.sum_depth = self.sum_depth + other.sum_depth
      result.min_depth = {self.min_depth, other.min_depth}.min
      result.max_depth = {self.max_depth, other.max_depth}.max
      result.bases = self.bases + other.bases
      result
    end
  end
end
