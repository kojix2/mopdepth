module Depth::Core
  module Cigar
    # Back-compatible helper that allocates a new array (prefer cigar_fill_events!)
    def cigar_start_end_events(cigar, ipos : Int32) : Array(Tuple(Int32, Int32))
      events = [] of Tuple(Int32, Int32)
      cigar_fill_events!(cigar, ipos, events)
      events
    end

    # Bit classification: bit0 = consumes reference, bit1 = consumes query
    private def classify_cigar(op : Char) : UInt8
      case op
      when 'M', '=', 'X' then 0b11_u8
      when 'D', 'N'      then 0b01_u8
      when 'I', 'S'      then 0b10_u8
      else                    0_u8
      end
    end

    # Streaming iterator for start/end events (pos, +1/-1) without allocations
    def cigar_each_event(cigar, ipos : Int32, & : Int32, Int32 ->)
      pos = ipos
      last_stop = -1
      cigar.each do |op_char, olen|
        cls = classify_cigar(op_char)
        next if (cls & 0b01) == 0
        len = olen.to_i32
        if (cls & 0b11) == 0b11
          if pos != last_stop
            yield pos, 1
            yield last_stop, -1 if last_stop >= 0
          end
          last_stop = pos + len
        end
        pos += len
      end
      yield last_stop, -1 if last_stop >= 0
    end

    # CIGAR → start/end events on reference, filling provided buffer
    def cigar_fill_events!(cigar, ipos : Int32, evbuf : Array(Tuple(Int32, Int32))) : Nil
      evbuf.clear
      cigar_each_event(cigar, ipos) do |p, v|
        evbuf << {p, v}
      end
      nil
    end

    def inc_coverage(cigar, ipos : Int32, a : Coverage)
      return if a.empty?
      cigar_each_event(cigar, ipos) do |pos, val|
        next if pos < 0
        p = pos.clamp(0, a.size - 1)
        a[p] += val
      end
    end

    # Build contiguous reference-aligned segments [start, stop) for M/= /X parts.
    def cigar_segments(cigar, ipos : Int32) : Array(Tuple(Int32, Int32))
      segs = [] of Tuple(Int32, Int32)
      pos = ipos
      seg_start = -1
      seg_stop = -1
      cigar.each do |op_char, olen|
        cls = classify_cigar(op_char)
        next if (cls & 0b01) == 0
        len = olen.to_i32
        if (cls & 0b11) == 0b11
          if pos != seg_stop
            segs << {seg_start, seg_stop} if seg_start >= 0 && seg_stop >= 0
            seg_start = pos
          end
          seg_stop = pos + len
        end
        pos += len
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
        cls = classify_cigar(op_char)
        next if (cls & 0b01) == 0
        len = olen.to_i32
        if (cls & 0b11) == 0b11
          if pos != seg_stop
            if seg_start >= 0 && seg_stop >= 0
              yield({seg_start, seg_stop})
            end
            seg_start = pos
          end
          seg_stop = pos + len
        end
        pos += len
      end
      if seg_start >= 0 && seg_stop >= 0
        yield({seg_start, seg_stop})
      end
    end
  end
end
