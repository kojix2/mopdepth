require "../core/region"
require "../core/interval"
require "../core/target"
require "../errors"
require "hts"
require "set"

module Depth::FileIO
  # Parse region like: chr1, chr1:100-200, or BED line
  def self.parse_region_str(s : String, targets : Array(Depth::Core::Target)) : Depth::Core::Region?
    return if s.empty?
    if target = targets.find { |candidate| candidate.name == s }
      return Depth::Core::Region.new(target.name, 0, 0)
    end
    return if s == "nil"

    target = targets.sort_by { |candidate| -candidate.name.bytesize }
      .find { |candidate| s.starts_with?(candidate.name + ":") }
    raise Depth::ConfigError.new("Chromosome not found or malformed region: #{s}") unless target

    coordinates = s[(target.name.bytesize + 1)..]
    match = /\A([0-9]+)(?:-([0-9]+))?\z/.match(coordinates)
    raise Depth::ConfigError.new("Malformed region: #{s}") unless match

    first = match[1].to_i64?
    last = (match[2]? || match[1]).to_i64?
    unless first && last && first >= 1 && last >= first && last <= target.length
      raise Depth::ConfigError.new("Region is outside #{target.name} (length #{target.length}): #{s}")
    end

    Depth::Core::Region.new(target.name, (first - 1).to_i32, last.to_i32)
  end

  # BED reader → {chrom => [Interval]}; intervals do not retain chromosome strings.
  def self.read_bed(path : String, targets : Array(Depth::Core::Target)) : Hash(String, Array(Depth::Core::Interval))
    tbl = Hash(String, Array(Depth::Core::Interval)).new { |hash, key| hash[key] = [] of Depth::Core::Interval }
    targets_by_name = targets.to_h { |target| {target.name, target} }
    unknown = Set(String).new
    line_number = 0

    HTS::Bgzf.open(path, "r") do |bed|
      bed.each_line do |line|
        line_number += 1
        next if line.strip.empty? || line.starts_with?("#") || line.starts_with?("track ") || line.starts_with?("browser ")
        if parsed = parse_bed_line(line, path, line_number, targets_by_name, unknown)
          chrom, interval = parsed
          tbl[chrom] << interval
        end
      end
    end

    tbl.each_value &.sort_by!(&.start)
    tbl
  end

  private def self.parse_bed_line(line : String, path : String, line_number : Int32,
                                  targets : Hash(String, Depth::Core::Target), unknown : Set(String)) : Tuple(String, Depth::Core::Interval)?
    cols = line.rstrip.split('\t', 5)
    raise Depth::ConfigError.new("#{path}:#{line_number}: BED line has fewer than 3 columns") if cols.size < 3

    chrom = cols[0]
    unless target = targets[chrom]?
      warn_unknown_chromosome(chrom, unknown)
      return
    end
    chrom_len = target.length

    start_pos = cols[1].to_i32?
    stop_pos = cols[2].to_i32?
    unless start_pos && stop_pos
      raise Depth::ConfigError.new("#{path}:#{line_number}: BED coordinates must be Int32 integers")
    end
    if start_pos < 0 || stop_pos <= start_pos || stop_pos > chrom_len
      raise Depth::ConfigError.new(
        "#{path}:#{line_number}: BED interval #{chrom}:#{start_pos}-#{stop_pos} is outside reference length #{chrom_len}"
      )
    end

    name = cols[3]?.try(&.strip).presence
    {target.name, Depth::Core::Interval.new(start_pos, stop_pos, name)}
  end

  private def self.warn_unknown_chromosome(chrom : String, unknown : Set(String))
    return if unknown.includes?(chrom)

    STDERR.puts "[mopdepth] warning: BED chromosome not found in alignment header: #{chrom}"
    unknown << chrom
  end
end
