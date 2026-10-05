require "option_parser"
require "./depth/config"
require "./depth/runner"
require "./depth/version"
require "./depth/errors"

module Depth
  class CLI
    def self.run(args = ARGV)
      Runner.new(parse_config(args)).run
    rescue ex : ConfigError | ArgumentError
      STDERR.puts "Config error: #{ex.message}"
      exit 1
    rescue ex : FileNotFoundError
      STDERR.puts "File not found: #{ex.message}"
      STDERR.puts ex.backtrace.first
      exit 1
    rescue ex : BamIndexError
      STDERR.puts "Index error: #{ex.message}"
      STDERR.puts ex.backtrace.first
      exit 1
    rescue ex : Exception
      STDERR.puts "Error: #{ex.message}"
      STDERR.puts ex.backtrace.join("\n")
      exit 1
    end

    private def self.parse_config(args : Array(String)) : Config
      config = Config.new
      argv = args.dup

      OptionParser.parse(argv) do |psr|
        psr.banner = "Usage: mopdepth [options] <prefix> <BAM-or-CRAM>"

        psr.on("-t", "--threads THREADS", "BAM decompression threads") { |v| config.threads = parse_i32(v, "--threads") }
        psr.on("-c", "--chrom CHROM", "Restrict to chromosome or 1-based range") { |v| config.chrom = v }
        psr.on("-b", "--by BY", "BED file or positive numeric window size") { |v| config.by = v }
        psr.on("-f", "--fasta FASTA", "FASTA reference for CRAM input [default: $REF_PATH]") { |v| config.fasta = v }
        psr.on("-n", "--no-per-base", "Skip per-base output") { config.no_per_base = true }
        psr.on("-Q", "--mapq MAPQ", "MAPQ threshold") { |v| config.mapq = parse_i32(v, "--mapq") }
        psr.on("-l", "--min-frag-len MIN", "Minimum fragment length") { |v| config.min_frag_len = parse_i32(v, "--min-frag-len") }
        psr.on("-u", "--max-frag-len MAX", "Maximum fragment length") { |v| config.max_frag_len = parse_i32(v, "--max-frag-len") }
        psr.on("-x", "--fast-mode", "Fast mode") { config.fast_mode = true }
        psr.on("-a", "--fragment-mode", "Count full fragment (proper pairs only)") { config.fragment_mode = true }
        psr.on("-m", "--use-median", "Use median for region stats instead of mean") { config.use_median = true }
        psr.on("-q", "--quantize QUANTIZE", "Write quantized output (e.g., 0:1:4:)") { |v| config.quantize = v }
        psr.on("-T", "--thresholds THRESHOLDS", "Comma-separated thresholds for region coverage") { |v| config.thresholds_str = v }
        psr.on("-F", "--flag FLAG", "Exclude reads with FLAG bits set") { |v| config.exclude_flag = parse_u16(v, "--flag") }
        psr.on("-i", "--include-flag FLAG", "Include only reads with FLAG bits set") { |v| config.include_flag = parse_u16(v, "--include-flag") }
        psr.on("-R", "--read-groups GROUPS", "Comma-separated read group IDs") { |v| config.read_groups_str = v }
        psr.on("-v", "--version", "Show version") { puts Depth::VERSION; exit 0 }
        psr.on("-M", "--mos", "Use mosdepth-compatible filenames (mosdepth.*)") { config.mos_style = true }
        psr.on("-h", "--help", "Show this message") { puts psr; exit 0 }
        psr.invalid_option { |opt| raise ConfigError.new("unknown option: #{opt}") }
      end

      unless argv.size == 2
        raise ConfigError.new("expected exactly <prefix> <BAM-or-CRAM>; got #{argv.size} positional arguments")
      end
      config.prefix = argv[0]
      config.path = argv[1]
      config
    end

    private def self.parse_i32(value : String, option : String) : Int32
      value.to_i32? || raise ConfigError.new("invalid integer for #{option}: #{value.inspect}")
    end

    private def self.parse_u16(value : String, option : String) : UInt16
      value.to_u16? || raise ConfigError.new("invalid integer for #{option}: #{value.inspect}")
    end
  end
end

Depth::CLI.run
