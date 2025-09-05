module Depth::Stats
  # Histogram-based median for integers
  class IntHistogram
    getter counts : Array(Int32)
    getter n : Int32 = 0
    # Track which bins were touched to avoid full-array zeroing on clear
    @mark : Array(Bool)
    @touched : Array(Int32)

    def initialize(size : Int32 = 65_536)
      @counts = size > 0 ? Array(Int32).new(size, 0) : [] of Int32
      @mark = size > 0 ? Array(Bool).new(size, false) : [] of Bool
      @touched = [] of Int32
      @n = 0
    end

    def add(value : Int32)
      @n += 1
      v = value < 0 ? raise ArgumentError.new("negative depth: #{value}") : value
      idx = v < counts.size ? v : counts.size - 1
      unless @mark[idx]
        @mark[idx] = true
        @touched << idx
      end
      counts[idx] += 1
    end

    def median : Int32
      return 0 if n == 0
      stop_n = ((n.to_f * 0.5) + 0.5).to_i
      cum = 0
      counts.each_with_index do |cnt, i|
        cum += cnt
        return i.to_i if cum >= stop_n
      end
      0
    end

    def clear
      return if n == 0
      @n = 0
      # zero only touched bins and reset marks
      @touched.each do |i|
        counts[i] = 0
        @mark[i] = false
      end
      @touched.clear
    end
  end
end
