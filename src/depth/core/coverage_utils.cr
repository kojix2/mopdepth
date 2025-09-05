module Depth::Core
  module CoverageUtils
    # Utility function for cumulative sum
    def prefix_sum!(a : Coverage)
      sum = 0
      i = 0
      while i < a.size
        sum += a[i]
        a[i] = sum
        i += 1
      end
    end

    # Yield intervals [start, stop) with constant depth value
    def each_constant_segment(a : Coverage, stop_at : Int32 = -1, &)
      return if a.size <= 1

      last_depth = Int32::MIN
      last_i = 0
      i = 0
      stop = stop_at <= 0 ? a.size - 1 : stop_at

      while i < stop
        depth = a[i]
        if depth == last_depth
          i += 1
          next
        end

        yield({last_i, i, last_depth}) if last_depth != Int32::MIN
        last_depth = depth
        last_i = i

        break if i + 1 == stop
        i += 1
      end

      # Final segment: cap at requested stop (exclusive)
      if last_i < stop
        yield({last_i, stop, last_depth})
      elsif last_i != i
        yield({last_i - 1, i, last_depth})
      else
        yield({last_i, i, last_depth})
      end
    end
  end
end
