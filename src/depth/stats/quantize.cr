module Depth::Stats
  module Quantize
    # Parse quantize arguments string into array of integers
    # Examples:
    #   ":1" -> [0, 1]
    #   "0:1:4:" -> [0, 1, 4, Int32::MAX]
    def self.get_quantize_args(qa : String) : Array(Int32)
      return [] of Int32 if qa == "nil" || qa.empty?

      a = qa

      # If no colons, wrap with colons
      if a.count(':') == 0
        a = ':' + a + ':'
      end

      # If starts with :, prepend 0
      if a[0] == ':'
        a = "0" + a
      end

      # If ends with :, append high value
      if a[-1] == ':'
        a = a + Int32::MAX.to_s
      end

      begin
        qs = a.split(':').map(&.to_i)
        qs.sort!
        qs
      rescue ex
        STDERR.puts "[mopdepth] invalid quantize string: '#{qa}'"
        exit(2)
      end
    end

    # Create lookup table for quantize bins
    # Examples:
    #   [0, 1, 4] -> ["0:1", "1:4", "4:inf"]
    def self.make_lookup(quants : Array(Int32)) : Array(String)
      return [] of String if quants.size <= 1

      lookup = [] of String

      (0...(quants.size - 1)).each do |i|
        # Check for custom environment variable labels
        env_var = "MOSDEPTH_Q#{i}"
        custom_label = ENV[env_var]?

        if custom_label
          lookup << custom_label
        else
          if quants[i + 1] == Int32::MAX
            lookup << "#{quants[i]}:inf"
          else
            lookup << "#{quants[i]}:#{quants[i + 1]}"
          end
        end
      end

      lookup
    end

    # Linear search to find which bin a value belongs to
    # Returns -1 if value is outside all bins
    def self.linear_search(q : Int32, vals : Array(Int32)) : Int32
      return -1 if vals.empty?
      return -1 if q < vals[0] || q > vals[-1]

      vals.each_with_index do |val, i|
        if val > q
          return i - 1
        elsif val == q
          return i
        end
      end

      vals.size - 1
    end

    # Generate quantized depth segments
    # limit: number of elements from coverage to consider (typically target_size)
    # Yields tuples of (start, stop, label) over [0, limit-1]
    def self.gen_quantized(quants : Array(Int32), coverage : Array(Int32), limit : Int32, & : Tuple(Int32, Int32, String) ->)
      return if quants.empty?
      return if coverage.empty?
      return if limit <= 1

      lookup = make_lookup(quants)
      return if lookup.empty?

      # Bound the iteration to the requested limit (effective length + 1)
      lim = Math.min(limit, coverage.size)

      last_quantized = linear_search(coverage[0], quants)
      last_pos = 0

      (0...(lim - 1)).each do |pos|
        depth = coverage[pos]
        quantized = linear_search(depth, quants)

        next if quantized == last_quantized

        if last_quantized != -1 && last_quantized < lookup.size
          yield({last_pos, pos, lookup[last_quantized]})
        end

        last_quantized = quantized
        last_pos = pos
      end

      # Handle the final segment
      if last_quantized != -1 && last_pos < lim - 1 && last_quantized < lookup.size
        yield({last_pos, lim - 1, lookup[last_quantized]})
      end
    end

    # Backward-compatible overload: consume entire coverage array
    def self.gen_quantized(quants : Array(Int32), coverage : Array(Int32), &block : Tuple(Int32, Int32, String) ->)
      gen_quantized(quants, coverage, coverage.size, &block)
    end
  end
end
