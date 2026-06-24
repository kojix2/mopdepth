require "hts"
require "./cigar"
require "./constants"
require "./coverage_result"
require "./coverage_utils"

module Depth::Core
  class CoverageCalculator
    include Cigar

    private class SeenMate
      getter start : Int32
      getter stop : Int32
      getter cigar_size : UInt32
      getter events : Array(Tuple(Int32, Int32))?

      def initialize(@start : Int32, @stop : Int32, @cigar_size : UInt32, @events : Array(Tuple(Int32, Int32))?)
      end
    end

    @evbuf : Array(Tuple(Int32, Int32))

    def initialize(@bam : HTS::Bam, @options : Options)
      # store first-read of overlapping proper pairs to correct double-counting when mate arrives
      @seen = Hash(String, SeenMate).new
      @evbuf = [] of Tuple(Int32, Int32)
    end

    # Check if record should be filtered out
    private def filtered_out?(rec) : Bool
      return true if rec.mapq < @options.mapq

      if @options.fragment_mode
        return true if @options.min_frag_len >= 0 && rec.isize.abs < @options.min_frag_len
        return true if rec.isize.abs > @options.max_frag_len
      end

      flag = rec.flag_value
      return true if (flag & @options.exclude_flag) != 0
      return true if @options.include_flag != 0 && (flag & @options.include_flag) == 0

      unless @options.read_groups.empty?
        rg = rec.aux_string("RG")
        return true unless rg && @options.read_groups.includes?(rg)
      end

      false
    end

    # Calculate coverage for a single record (positions are shifted by `offset` for region queries)

    private def accumulate_record!(rec, coverage : Coverage, offset : Int32)
      if @options.fast_mode
        start_pos = (rec.pos.to_i32 - offset).clamp(0, coverage.size - 1)
        mark_and_add!(coverage, start_pos, 1)
        endp = (rec.endpos.to_i32 - offset)
        endp = coverage.size - 1 if endp >= coverage.size
        endp = 0 if endp < 0
        mark_and_add!(coverage, endp, -1)
      elsif @options.fragment_mode
        return if rec.read2? || !rec.proper_pair? || rec.supplementary?
        frag_start = Math.min(rec.pos, rec.mate_pos).to_i32 - offset
        frag_len = rec.isize.abs
        end_pos = frag_start + frag_len
        end_pos = coverage.size - 1 if end_pos >= coverage.size
        startp = frag_start.clamp(0, coverage.size - 1)
        mark_and_add!(coverage, startp, 1)
        mark_and_add!(coverage, end_pos, -1)
      else
        # Default (per-base) mode with mosdepth-like mate-overlap correction
        if @options.fast_mode == false && @options.fragment_mode == false &&
           rec.proper_pair? && !rec.supplementary? && rec.tid == rec.mtid && rec.mate_pos >= 0
          rec_start = rec.pos.to_i32
          rec_stop = rec.endpos.to_i32
          # If this read overlaps its mate and is the earlier (or equal) one, store it; otherwise, if mate was stored, correct overlap now
          if rec_stop > rec.mate_pos
            qname = rec.qname
            if rec_start < rec.mate_pos || (rec_start == rec.mate_pos && !@seen.has_key?(qname))
              @seen[qname] = capture_seen_mate(rec, rec_start, rec_stop)
            elsif mate = @seen.delete(qname)
              correct_mate_overlap_coverage!(rec, mate, rec_start, rec_stop, coverage, offset)
            end
          else
            if mate = @seen.delete(rec.qname)
              correct_mate_overlap_coverage!(rec, mate, rec_start, rec_stop, coverage, offset)
            end
          end
        end
        # Always add coverage for this record after possible overlap correction step
        apply_record_cigar!(rec, rec.pos.to_i32, coverage, offset)
      end
    end

    private def capture_seen_mate(rec, rec_start : Int32, rec_stop : Int32) : SeenMate
      cigar_size = rec.cigar_size
      events = nil
      if cigar_size != 1
        ev = [] of Tuple(Int32, Int32)
        record_cigar_append_events!(rec, rec_start, ev)
        events = ev
      end
      SeenMate.new(rec_start, rec_stop, cigar_size, events)
    end

    private def append_seen_mate_events!(mate : SeenMate, events : Array(Tuple(Int32, Int32)))
      if mate_events = mate.events
        events.concat(mate_events)
      else
        events << {mate.start, 1}
        events << {mate.stop, -1}
      end
    end

    private def correct_mate_overlap_coverage!(rec, mate : SeenMate, rec_start : Int32, rec_stop : Int32,
                                               coverage : Coverage, offset : Int32)
      if rec.cigar_size == 1 && mate.cigar_size == 1
        s = [rec_start, mate.start].max
        e = [rec_stop, mate.stop].min
        if e > s
          s = (s - offset).clamp(0, coverage.size - 1)
          e = (e - offset).clamp(0, coverage.size - 1)
          mark_and_add!(coverage, s, -1)
          mark_and_add!(coverage, e, 1)
        end
      else
        ses = @evbuf
        ses.clear
        record_cigar_append_events!(rec, rec_start, ses)
        append_seen_mate_events!(mate, ses)
        ses.sort! { |a, b| a[0] <=> b[0] }
        pair_depth = 0
        last_pos = 0
        ses.each do |pos, val|
          if val == -1 && pair_depth == 2
            s = last_pos
            e = pos
            if e > s
              s = (s - offset).clamp(0, coverage.size - 1)
              e = (e - offset).clamp(0, coverage.size - 1)
              mark_and_add!(coverage, s, -1)
              mark_and_add!(coverage, e, 1)
            end
          end
          pair_depth += val
          last_pos = pos
        end
      end
    end

    # Apply cigar events into diff-array
    private def apply_cigar!(cigar, ipos : Int32, a : Coverage, offset : Int32)
      @evbuf.clear
      cigar_fill_events!(cigar, ipos, @evbuf)
      @evbuf.each do |pos, val|
        p = (pos - offset).clamp(0, a.size - 1)
        mark_and_add!(a, p, val)
      end
    end

    private def apply_record_cigar!(rec, ipos : Int32, a : Coverage, offset : Int32)
      record_cigar_each_event(rec, ipos) do |pos, val|
        p = (pos - offset).clamp(0, a.size - 1)
        mark_and_add!(a, p, val)
      end
    end

    # Initialize coverage array for a chromosome
    def initialize_coverage_array(coverage : Coverage, chrom_len : Int32)
      target_size = chrom_len + 1
      if coverage.size == target_size
        # caller should reset; avoid full memset here
      else
        coverage.clear
        coverage.concat(Array(Int32).new(target_size, 0))
      end
    end

    # Process region-specific query
    private def process_region_query(a : Coverage, r : Region, tid : Int32, offset : Int32) : {Bool, Int32}
      found = false
      chrom_tid = UNKNOWN_CHROM_TID

      q_start = (r.start > 0 ? r.start : 0)
      q_stop = (r.stop > 0 ? r.stop : @bam.header.target_len[tid].to_i32)

      # Assume Runner already sized/zeroed the coverage buffer; avoid duplicate initialization here
      @seen.clear
      @bam.query(tid, q_start, q_stop) do |rec|
        next if filtered_out?(rec)
        found = true unless found
        accumulate_record!(rec, a, offset)
      end

      {found, chrom_tid}
    end

    # Process full BAM scan
    private def process_full_scan(a : Coverage, tid : Int32) : {Bool, Int32}
      found = false
      chrom_tid = UNKNOWN_CHROM_TID
      current_tid = Int32::MIN

      @seen.clear
      @bam.each(copy: false) do |rec|
        # Break if we encounter a different chromosome to prevent array corruption
        if current_tid != Int32::MIN && rec.tid != current_tid
          break
        end

        unless found
          chrom_tid = (tid >= 0) ? tid : rec.tid
          current_tid = chrom_tid
          found = true
        end

        next if filtered_out?(rec)
        accumulate_record!(rec, a, 0)
      end

      {found, chrom_tid}
    end

    # Returns: CoverageResult enum values or tid
    def calculate(a : Coverage, r : Region?, offset : Int32 = 0) : Int32
      # Determine tid
      tid = CoverageResult::ChromNotFound.value
      if r
        tid = @bam.header.get_tid(r.chrom)
      else
        # if no region chrom provided, process everything by iteration
        tid = CoverageResult::AllChroms.value
      end

      return CoverageResult::ChromNotFound.value if tid == CoverageResult::ChromNotFound.value

      found, chrom_tid = if r && tid >= 0
                           process_region_query(a, r, tid, offset)
                         else
                           process_full_scan(a, tid)
                         end

      return CoverageResult::NoData.value unless found
      tid >= 0 ? tid : chrom_tid
    end

    private def mark_and_add!(a : Coverage, idx : Int32, delta : Int32)
      a[idx] += delta
    end
  end
end
