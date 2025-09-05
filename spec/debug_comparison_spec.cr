require "./spec_helper"
require "file_utils"
require "process"

def mosdepth_available_for_debug?(mosdepth_dir : String) : Bool
  mosdepth_path = ENV.fetch("MOSDEPTH_PATH", File.expand_path("#{mosdepth_dir}/mosdepth"))
  begin
    status = Process.run(mosdepth_path, ["-h"],
      output: Process::Redirect::Close,
      error: Process::Redirect::Close)
    status.success?
  rescue
    false
  end
end

def run_mosdepth_debug(args : Array(String), prefix : String, temp_dir : String, mosdepth_dir : String, test_bam : String)
  mosdepth_path = ENV.fetch("MOSDEPTH_PATH", File.expand_path("#{mosdepth_dir}/mosdepth"))
  test_bam_path = File.expand_path(test_bam)
  full_args = [prefix] + args + [test_bam_path]

  puts "Running mosdepth: #{mosdepth_path} #{full_args.join(" ")}"
  puts "Working directory: #{temp_dir}"
  puts "BAM file exists: #{File.exists?(test_bam_path)}"
  puts "mosdepth binary exists: #{File.exists?(mosdepth_path)}"

  stdout = IO::Memory.new
  stderr = IO::Memory.new

  result = Process.run(mosdepth_path, full_args,
    chdir: temp_dir,
    output: stdout,
    error: stderr)

  puts "Exit code: #{result.exit_code}"
  puts "STDOUT: #{stdout}"
  puts "STDERR: #{stderr}"
  puts "---"

  result
end

def run_mopdepth_debug(args : Array(String), prefix : String, temp_dir : String, test_bam : String)
  full_args = args + [prefix, test_bam]

  bin = TestBin.binary
  puts "Running mopdepth: #{bin} #{(args + [prefix, test_bam]).join(" ")}"
  puts "Working directory: #{Dir.current}"
  puts "Output directory: #{temp_dir}"

  stdout = IO::Memory.new
  stderr = IO::Memory.new

  TestBin.ensure_built!
  result = Process.run(bin, args + [prefix, test_bam],
    chdir: Dir.current,
    env: {"PWD" => temp_dir},
    output: stdout,
    error: stderr)

  puts "Exit code: #{result.exit_code}"
  puts "STDOUT: #{stdout}"
  puts "STDERR: #{stderr}"
  puts "---"

  result
end

describe "Debug comparison with mosdepth" do
  temp_dir = "/tmp/mopdepth_debug"
  mosdepth_dir = "./mosdepth"
  test_bam = "#{mosdepth_dir}/tests/ovl.bam"

  before_each do
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "debugs basic coverage calculation" do
    # Skip this diagnostic test if mosdepth binary is not available
    next unless mosdepth_available_for_debug?(mosdepth_dir)
    puts "\n=== DEBUGGING BASIC COVERAGE CALCULATION ==="

    puts "\n--- Running mosdepth ---"
    mosdepth_result = run_mosdepth_debug([] of String, "#{temp_dir}/mosdepth_test", temp_dir, mosdepth_dir, test_bam)

    puts "\n--- Running mopdepth ---"
    mopdepth_result = run_mopdepth_debug([] of String, "#{temp_dir}/mopdepth_test", temp_dir, test_bam)

    puts "\n--- File listing after execution ---"
    if Dir.exists?(temp_dir)
      Dir.each_child(temp_dir) do |file|
        full_path = File.join(temp_dir, file)
        size = File.exists?(full_path) ? File.size(full_path) : 0
        puts "#{file}: #{size} bytes"
      end
    end

    # Don't fail the test, just show the results
    puts "\nmosdepth success: #{mosdepth_result.success?}"
    puts "mopdepth success: #{mopdepth_result.success?}"
  end
end
