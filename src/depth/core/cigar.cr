module Depth::Core
  module Cigar
    # Back-compatible helper that allocates a new array (prefer cigar_fill_events!)
    def cigar_start_end_events(cigar, ipos : Int32) : Array(Tuple(Int32, Int32))
      events = [] of Tuple(Int32, Int32)
      cigar_fill_events!(cigar, ipos, events)
      events
    end

    # CIGAR → start/end events on reference
    # Returns array of {pos, +1/-1}; pos is 0-based ref coordinate
    #
    # Key points (what contributes to depth):
    # - Only operations that consume both reference and query contribute to depth: M, =, X
    # - D, N consume reference only (gaps). They do not add to depth; we just advance the ref position (pos)
    # - I, S consume query only. They do not advance reference and do not add to depth
    #
    # Implementation strategy:
    # - Operations that don't consume reference are skipped entirely (next unless consumes_ref)
    # - Only when an op consumes both reference and query (M/= /X) we emit start(+1)/end(-1) events
    #   This ensures only M/= /X spans change coverage when later applied to the diff array
    # Fill start/end events into a provided buffer to avoid allocations.
    # Each event is a tuple: {pos, +1/-1}
    def cigar_fill_events!(cigar, ipos : Int32, evbuf : Array(Tuple(Int32, Int32))) : Nil
      pos = ipos
      last_stop = -1
      cigar.each do |op_char, olen|
        # Check if operation consumes reference: M, D, N, =, X
        consumes_ref = ['M', 'D', 'N', '=', 'X'].includes?(op_char)
        next unless consumes_ref

        # Check if operation consumes query: M, I, S, =, X
        consumes_query = ['M', 'I', 'S', '=', 'X'].includes?(op_char)
        if consumes_query
          # We get here only for M/= /X (I/S don't consume reference and were filtered above).
          # If the previous segment doesn't continue, start a new one;
          # if a prior segment exists, emit its end event.
          if pos != last_stop
            evbuf << {pos, 1}
            evbuf << {last_stop, -1} if last_stop >= 0
          end
          last_stop = pos + olen.to_i32
        end
        # For D/N (consume reference only), consumes_query is false:
        # don't create events, only advance pos (treat as a gap)
        pos += olen.to_i32
      end
      # Close the last open segment if any
      evbuf << {last_stop, -1} if last_stop >= 0
      nil
    end

    def inc_coverage(cigar, ipos : Int32, a : Coverage)
      # Apply the start(+1)/end(-1) events to the diff array 'a'.
      # A subsequent prefix_sum! turns it into actual per-base coverage.
      cigar_start_end_events(cigar, ipos).each do |pos, val|
        next if pos < 0 || a.empty?
        p = pos.clamp(0, a.size - 1)
        a[p] += val
      end
    end

    # Build contiguous reference-aligned segments [start, stop) for M/= /X parts.
    # This mirrors cigar_start_end_events but returns merged intervals directly.
    def cigar_segments(cigar, ipos : Int32) : Array(Tuple(Int32, Int32))
      segs = [] of Tuple(Int32, Int32)
      pos = ipos
      seg_start = -1
      seg_stop = -1
      cigar.each do |op_char, olen|
        consumes_ref = ['M', 'D', 'N', '=', 'X'].includes?(op_char)
        next unless consumes_ref

        consumes_query = ['M', 'I', 'S', '=', 'X'].includes?(op_char)
        if consumes_query
          # start new segment if discontinuous
          if pos != seg_stop
            # close previous
            segs << {seg_start, seg_stop} if seg_start >= 0 && seg_stop >= 0
            seg_start = pos
          end
          seg_stop = pos + olen.to_i32
        end
        pos += olen.to_i32
      end
      segs << {seg_start, seg_stop} if seg_start >= 0 && seg_stop >= 0
      segs
    end

    # Stream segments without building an array
    def cigar_each_segment(cigar, ipos : Int32, & : Tuple(Int32, Int32) ->)
      pos = ipos
      seg_start = -1
      seg_stop = -1
      cigar.each do |op_char, olen|
        consumes_ref = ['M', 'D', 'N', '=', 'X'].includes?(op_char)
        unless consumes_ref
          next
        end
        consumes_query = ['M', 'I', 'S', '=', 'X'].includes?(op_char)
        if consumes_query
          if pos != seg_stop
            if seg_start >= 0 && seg_stop >= 0
              yield({seg_start, seg_stop})
            end
            seg_start = pos
          end
          seg_stop = pos + olen.to_i32
        end
        pos += olen.to_i32
      end
      if seg_start >= 0 && seg_stop >= 0
        yield({seg_start, seg_stop})
      end
    end
  end
end
