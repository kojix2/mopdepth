require "hts"
require "./cigar"
require "./constants"
require "./coverage_result"
require "./coverage_utils"

module Depth::Core
  class CoverageCalculator
    include Cigar
    # touched tracking
    @generation : UInt32
    @marks : Array(UInt32)
    @touched : Array(Int32)
    @evbuf : Array(Tuple(Int32, Int32))

    def initialize(@bam : HTS::Bam, @options : Options)
      # store first-read of overlapping proper pairs to correct double-counting when mate arrives
      @seen = Hash(String, HTS::Bam::Record).new
      # generation-based touched tracking (avoid full memset)
      @generation = 1_u32
      @marks = [] of UInt32  # same length as coverage capacity
      @touched = [] of Int32 # indices touched during current generation
      @evbuf = [] of Tuple(Int32, Int32)
    end

    # Check if record should be filtered out
    private def filtered_out?(rec) : Bool
      return true if rec.mapq < @options.mapq

      if @options.fragment_mode
        return true if @options.min_frag_len >= 0 && rec.isize.abs < @options.min_frag_len
        return true if rec.isize.abs > @options.max_frag_len
      end

      return true if (rec.flag.value & @options.exclude_flag) != 0
      return true if @options.include_flag != 0 && (rec.flag.value & @options.include_flag) == 0

      unless @options.read_groups.empty?
        rg = rec.aux("RG")
        return true unless rg.is_a?(String) && @options.read_groups.includes?(rg)
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
        return if rec.flag.read2? || !rec.flag.proper_pair? || rec.flag.supplementary?
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
           rec.flag.proper_pair? && !rec.flag.supplementary? && rec.tid == rec.mtid && rec.mate_pos >= 0
          rec_start = rec.pos.to_i32
          rec_stop = rec.endpos.to_i32
          # If this read overlaps its mate and is the earlier (or equal) one, store it; otherwise, if mate was stored, correct overlap now
          if rec_stop > rec.mate_pos && (rec_start < rec.mate_pos || (rec_start == rec.mate_pos && !@seen.has_key?(rec.qname)))
            # store a clone since records are reused
            @seen[rec.qname] = rec.clone
          else
            if mate = @seen.delete(rec.qname)
              # Fast path: both reads have single M op
              if rec.cigar.size == 1 && mate.cigar.size == 1
                s = [rec_start, mate.pos.to_i32].max
                e = [rec_stop, mate.endpos.to_i32].min
                if e > s
                  s = (s - offset).clamp(0, coverage.size - 1)
                  e = (e - offset).clamp(0, coverage.size - 1)
                  mark_and_add!(coverage, s, -1)
                  mark_and_add!(coverage, e, 1)
                end
              else
                # Build combined start/end events for rec and mate, subtract where pair_depth==2
                ses = @evbuf
                ses.clear
                cigar_fill_events!(rec.cigar, rec_start, ses)
                cigar_fill_events!(mate.cigar, mate.pos.to_i32, ses)
                ses.sort_by! { |(p, _)| p }
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
          end
        end
        # Always add coverage for this record after possible overlap correction step
        apply_cigar!(rec.cigar, rec.pos.to_i32, coverage, offset)
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

    # Initialize coverage array for a chromosome
    def initialize_coverage_array(coverage : Coverage, chrom_len : Int32)
      target_size = chrom_len + 1
      if coverage.size == target_size
        # caller should reset; avoid full memset here
      else
        coverage.clear
        coverage.concat(Array(Int32).new(target_size, 0))
        ensure_marks_capacity(target_size)
      end
    end

    # Process region-specific query
    private def process_region_query(a : Coverage, r : Region, tid : Int32, offset : Int32) : {Bool, Int32}
      found = false
      chrom_tid = UNKNOWN_CHROM_TID

      q_start = (r.start > 0 ? r.start : 0)
      q_stop = (r.stop > 0 ? r.stop : @bam.header.target_len[tid].to_i32)
      region_str = "#{r.chrom}:#{q_start + 1}-#{q_stop}"

      # Assume Runner already sized/zeroed the coverage buffer; avoid duplicate initialization here
      @seen.clear
      @bam.query(region_str) do |rec|
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

    # Reset using touch marks: zero only indices written in current generation (up to limit)
    def reset_coverage!(a : Coverage, limit : Int32)
      if @touched.size > 0
        lim = Math.min(limit, a.size)
        @touched.each do |idx|
          next if idx >= lim
          a[idx] = 0
        end
        @touched.clear
      end
      # advance generation (lazy clear of marks)
      @generation &+= 1_u32
      if @generation == 0_u32
        # wrapped; hard reset marks
        @marks.fill(0_u32)
        @generation = 1_u32
      end
    end

    private def ensure_marks_capacity(capacity : Int32)
      if @marks.size < capacity
        @marks.concat(Array(UInt32).new(capacity - @marks.size, 0_u32))
      end
    end

    private def mark_and_add!(a : Coverage, idx : Int32, delta : Int32)
      ensure_marks_capacity(a.size)
      # mark once per generation
      if @marks[idx] != @generation
        @marks[idx] = @generation
        @touched << idx
        a[idx] = 0 if a[idx] != 0
      end
      a[idx] += delta
    end
  end
end
