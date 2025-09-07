require "./spec_helper"
require "../src/depth/core/cigar" # load Depth::Core::Cigar

# Compatibility tests ensuring optimized cigar code preserves behavior.
# These tests rely on public APIs: cigar_fill_events!, cigar_segments, cigar_each_segment.
# We construct synthetic CIGAR objects using the same interface as production code (each yielding (Char, Int)).

# Minimal helper fake CIGAR container for testing if real parser is not available here.
class TestCigar
  getter ops : Array(Tuple(Char, Int32))

  def initialize(@ops : Array(Tuple(Char, Int32))); end

  def each(&block : Char, Int32 ->)
    @ops.each { |op| yield op[0], op[1] }
  end
end

# Reference implementation (pre-optimization) for events and segments for comparison.
module ReferenceCigarImpl
  extend self

  def fill_events!(cigar, ipos : Int32, evbuf : Array(Tuple(Int32, Int32)))
    pos = ipos
    last_stop = -1
    cigar.each do |op_char, olen|
      consumes_ref = case op_char
                     when 'M', 'D', 'N', '=', 'X' then true
                     else                              false
                     end
      next unless consumes_ref
      consumes_query = case op_char
                       when 'M', 'I', 'S', '=', 'X' then true
                       else                              false
                       end
      if consumes_query
        if pos != last_stop
          evbuf << {pos, 1}
          evbuf << {last_stop, -1} if last_stop >= 0
        end
        last_stop = pos + olen.to_i32
      end
      pos += olen.to_i32
    end
    evbuf << {last_stop, -1} if last_stop >= 0
  end

  def segments(cigar, ipos : Int32)
    segs = [] of Tuple(Int32, Int32)
    pos = ipos
    seg_start = -1
    seg_stop = -1
    cigar.each do |op_char, olen|
      consumes_ref = case op_char
                     when 'M', 'D', 'N', '=', 'X' then true
                     else                              false
                     end
      next unless consumes_ref
      consumes_query = case op_char
                       when 'M', 'I', 'S', '=', 'X' then true
                       else                              false
                       end
      if consumes_query
        if pos != seg_stop
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
end

# Helper class including the optimized module under test.
class CigarHost
  include Depth::Core::Cigar
end

HOST = CigarHost.new

# helper outside describe to avoid dynamic def issue
private def build_test_cigar(str : String)
  ops = [] of Tuple(Char, Int32)
  num = 0
  str.each_char do |c|
    if c.ascii_number?
      num = num * 10 + (c.ord - 48)
    else
      raise "length missing before op #{c}" if num == 0
      ops << {c, num.to_i}
      num = 0
    end
  end
  TestCigar.new(ops)
end

describe "Cigar optimization compatibility" do
  it "events match reference" do
    cases = [
      "10M",
      "5M3I7M2D4M",
      "10M5S8M",
      "3S5M2N6M1D5M4I9M",
      "6M1I1M1D2M1N3M",
      "1M1D1M1N1M1I1M1S1M",
    ]
    cases.each do |c|
      cigar = build_test_cigar(c)
      (0..3).each do |ipos|
        ref_buf = [] of Tuple(Int32, Int32)
        new_buf = [] of Tuple(Int32, Int32)
        ReferenceCigarImpl.fill_events!(cigar, ipos, ref_buf)
        HOST.cigar_fill_events!(cigar, ipos, new_buf)
        ref_buf.should eq new_buf
      end
    end
  end

  it "segments match reference" do
    cases = [
      "10M",
      "5M3I7M2D4M",
      "10M5S8M",
      "3S5M2N6M1D5M4I9M",
      "6M1I1M1D2M1N3M",
      "1M1D1M1N1M1I1M1S1M",
    ]
    cases.each do |c|
      cigar = build_test_cigar(c)
      (0..3).each do |ipos|
        ref_segs = ReferenceCigarImpl.segments(cigar, ipos)
        new_segs = HOST.cigar_segments(cigar, ipos)
        ref_segs.should eq new_segs
      end
    end
  end

  it "cigar_each_event yields same sequence as fill_events!" do
    cigar = build_test_cigar("5M2D3M4M1N2M")
    buf = [] of Tuple(Int32, Int32)
    HOST.cigar_fill_events!(cigar, 0, buf)
    yielded = [] of Tuple(Int32, Int32)
    HOST.cigar_each_event(cigar, 0) do |p, v|
      yielded << {p, v}
    end
    yielded.should eq buf
  end
end
