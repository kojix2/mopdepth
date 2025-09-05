require "option_parser"
require "./depth/config"
require "./depth/runner"
require "./depth/version"
require "./depth/errors"

module Depth
  class CLI
    def self.run(args = ARGV)
      config = Config.new
      config.prefix = "out"

      OptionParser.parse(args) do |psr|
        psr.banner = "Usage: mopdepth [options] <prefix> <BAM-or-CRAM>"

        psr.on("-t", "--threads THREADS", "BAM decompression threads") { |v| config.threads = v.to_i }
        psr.on("-c", "--chrom CHROM", "Restrict to chromosome") { |v| config.chrom = v }
        psr.on("-b", "--by BY", "BED file or numeric window size") { |v| config.by = v }
        psr.on("-n", "--no-per-base", "Skip per-base output") { config.no_per_base = true }
        psr.on("-Q", "--mapq MAPQ", "MAPQ threshold") { |v| config.mapq = v.to_i }
        psr.on("-l", "--min-frag-len MIN", "Minimum fragment length") { |v| config.min_frag_len = v.to_i }
        psr.on("-u", "--max-frag-len MAX", "Maximum fragment length") { |v| config.max_frag_len = v.to_i }
        psr.on("-x", "--fast-mode", "Fast mode") { config.fast_mode = true }
        psr.on("-a", "--fragment-mode", "Count full fragment (proper pairs only)") { config.fragment_mode = true }
        psr.on("-m", "--use-median", "Use median for region stats instead of mean") { config.use_median = true }
        psr.on("-q", "--quantize QUANTIZE", "Write quantized output (e.g., 0:1:4:)") { |v| config.quantize = v }
        psr.on("-T", "--thresholds THRESHOLDS", "Comma-separated thresholds for region coverage") { |v| config.thresholds_str = v }
        psr.on("-F", "--flag FLAG", "Exclude reads with FLAG bits set") { |v| config.exclude_flag = v.to_u16 }
        psr.on("-i", "--include-flag FLAG", "Include only reads with FLAG bits set") { |v| config.include_flag = v.to_u16 }
        psr.on("-R", "--read-groups GROUPS", "Comma-separated read group IDs") { |v| config.read_groups_str = v }
        psr.on("-v", "--version", "Show version") { puts Depth::VERSION; exit 0 }
        psr.on("-M", "--mos", "Use mosdepth-compatible filenames (mosdepth.*)") { config.mos_style = true }
        psr.on("-h", "--help", "Show this message") { puts psr; exit 0 }

        psr.invalid_option do |opt|
          STDERR.puts "Error: unknown option: #{opt}"
          STDERR.puts "Use --help for usage information"
          exit 1
        end
      end

      if ARGV.size < 2
        STDERR.puts "Error: missing arguments: <prefix> <BAM-or-CRAM>"
        STDERR.puts "Use --help for usage information"
        exit 2
      end

      config.prefix = ARGV[0]
      config.path = ARGV[1]

      begin
        runner = Runner.new(config)
        runner.run
      rescue ex : ConfigError
        STDERR.puts "Config error: #{ex.message}"
        exit 1
      rescue ex : FileNotFoundError
        STDERR.puts "File not found: #{ex.message}"
        exit 1
      rescue ex : BamIndexError
        STDERR.puts "Index error: #{ex.message}"
        exit 1
      rescue ex : Exception
        STDERR.puts "Error: #{ex.message}"
        exit 1
      end
    end
  end
end

Depth::CLI.run
