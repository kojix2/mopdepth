require "../core/region"

module Depth::FileIO
  # Parse region like: chr1, chr1:100-200, or BED line
  def self.parse_region_str(s : String) : Depth::Core::Region?
    return nil if s.empty? || s == "nil"
    if s.includes?(':')
      chrom, rest = s.split(":", 2)
      if rest.includes?('-')
        st, ed = rest.split('-', 2)
        Depth::Core::Region.new(chrom, st.to_i - 1, ed.to_i)
      else
        Depth::Core::Region.new(chrom, rest.to_i - 1, rest.to_i)
      end
    else
      Depth::Core::Region.new(s, 0, 0)
    end
  end

  # BED reader → {chrom => [Region]}
  def self.read_bed(path : String) : Hash(String, Array(Depth::Core::Region))
    tbl = Hash(String, Array(Depth::Core::Region)).new { |h, k| h[k] = [] of Depth::Core::Region }
    File.each_line(path) do |line|
      next if line.starts_with?("#") || line.empty?
      next if line.starts_with?("track ")
      cols = line.rstrip.split('\t')
      if cols.size < 3
        STDERR.puts "[mopdepth] skipping bad bed line: #{line}"
        next
      end
      chrom = cols[0]
      s = cols[1].to_i
      e = cols[2].to_i
      name = cols.size > 3 ? cols[3].strip : nil
      tbl[chrom] << Depth::Core::Region.new(chrom, s, e, name)
    end
    # sort in-place
    tbl.each_value &.sort_by!(&.start)
    tbl
  end
end
