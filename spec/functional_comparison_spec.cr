require "./spec_helper"
require "file_utils"
require "process"

def mosdepth_available?(mosdepth_dir : String) : Bool
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

def run_mosdepth(args : Array(String), prefix : String, temp_dir : String, mosdepth_dir : String, test_bam : String) : Process::Status
  # Use MOSDEPTH_PATH environment variable if set, otherwise use local mosdepth
  mosdepth_path = ENV.fetch("MOSDEPTH_PATH", File.expand_path("#{mosdepth_dir}/mosdepth"))
  test_bam_path = File.expand_path(test_bam)
  full_args = [prefix] + args + [test_bam_path]

  Process.run(mosdepth_path, full_args,
    chdir: temp_dir,
    output: Process::Redirect::Close,
    error: Process::Redirect::Close)
end

def run_mopdepth(args : Array(String), prefix : String, temp_dir : String, test_bam : String) : Process::Status
  TestBin.run(["-M"] + args, prefix, test_bam, temp_dir)
end

def compare_files(file1 : String, file2 : String) : Bool
  return false unless File.exists?(file1) && File.exists?(file2)

  # Read content with per-file gzip handling
  content1 = TestIO.read_text(file1)
  content2 = TestIO.read_text(file2)

  content1.strip == content2.strip
end

describe "Functional comparison with mosdepth" do
  temp_dir = "/tmp/mopdepth_test"
  mosdepth_dir = "./mosdepth"
  test_bam = "#{mosdepth_dir}/tests/ovl.bam"

  before_each do
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  describe "basic coverage calculation" do
    it "produces identical per-base output" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth([] of String, "#{temp_dir}/mosdepth_test", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth([] of String, "#{temp_dir}/mopdepth_test", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      # Compare per-base output
      mosdepth_file = "#{temp_dir}/mosdepth_test.per-base.bed.gz"
      mopdepth_file = "#{temp_dir}/mopdepth_test.per-base.bed.gz"

      File.exists?(mosdepth_file).should be_true
      File.exists?(mopdepth_file).should be_true

      # Compare contents (handling gzip for mosdepth)
      compare_files(mosdepth_file, mopdepth_file).should be_true
    end

    it "produces identical summary output" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth([] of String, "#{temp_dir}/mosdepth_summary", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth([] of String, "#{temp_dir}/mopdepth_summary", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mosdepth_summary = "#{temp_dir}/mosdepth_summary.mosdepth.summary.txt"
      mopdepth_summary = "#{temp_dir}/mopdepth_summary.mosdepth.summary.txt"

      File.exists?(mosdepth_summary).should be_true
      File.exists?(mopdepth_summary).should be_true

      # Compare summary file contents
      compare_files(mosdepth_summary, mopdepth_summary).should be_true
    end
  end

  describe "flag filtering" do
    it "handles exclude flag filtering" do
      next unless mosdepth_available?(mosdepth_dir)
      # Test with -F 4 (exclude unmapped reads)
      mosdepth_status = run_mosdepth(["-F", "4"], "#{temp_dir}/mosdepth_flag", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-F", "4"], "#{temp_dir}/mopdepth_flag", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      # Both should produce output files
      mos_file = "#{temp_dir}/mosdepth_flag.per-base.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_flag.per-base.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end

    it "handles include flag filtering" do
      next unless mosdepth_available?(mosdepth_dir)
      # Test with -i 2 (include only proper pairs)
      mosdepth_status = run_mosdepth(["-i", "2"], "#{temp_dir}/mosdepth_include", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-i", "2"], "#{temp_dir}/mopdepth_include", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_include.per-base.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_include.per-base.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "quantization" do
    it "handles quantization identically" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["-q", "0:1:1000"], "#{temp_dir}/mosdepth_quant", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-q", "0:1:1000"], "#{temp_dir}/mopdepth_quant", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_quant.quantized.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_quant.quantized.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "region processing" do
    it "handles window-based regions" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["-b", "100"], "#{temp_dir}/mosdepth_window", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-b", "100"], "#{temp_dir}/mopdepth_window", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_window.regions.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_window.regions.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end

    it "handles BED file regions" do
      next unless mosdepth_available?(mosdepth_dir)
      bed_file = File.expand_path("#{mosdepth_dir}/tests/track.bed")
      next unless File.exists?(bed_file)

      mosdepth_status = run_mosdepth(["-b", bed_file], "#{temp_dir}/mosdepth_bed", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-b", bed_file], "#{temp_dir}/mopdepth_bed", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_bed.regions.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_bed.regions.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "threshold processing" do
    it "handles threshold calculations" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["-T", "0,1,2", "-b", "100"], "#{temp_dir}/mosdepth_thresh", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-T", "0,1,2", "-b", "100"], "#{temp_dir}/mopdepth_thresh", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_thresh.thresholds.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_thresh.thresholds.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "fast mode" do
    it "handles fast mode processing" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["-x"], "#{temp_dir}/mosdepth_fast", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-x"], "#{temp_dir}/mopdepth_fast", temp_dir, test_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_fast.per-base.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_fast.per-base.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "fragment mode" do
    it "handles fragment mode processing" do
      next unless mosdepth_available?(mosdepth_dir)
      fragment_bam = File.expand_path("#{mosdepth_dir}/tests/full-fragment-pairs.bam")
      next unless File.exists?(fragment_bam)

      # Use shared helpers to ensure MOSDEPTH_PATH is respected and paths are absolute
      mosdepth_status = run_mosdepth(["-a"], "#{temp_dir}/mosdepth_frag", temp_dir, mosdepth_dir, fragment_bam)
      mopdepth_status = run_mopdepth(["-a"], "#{temp_dir}/mopdepth_frag", temp_dir, fragment_bam)

      mosdepth_status.success?.should be_true
      mopdepth_status.success?.should be_true

      mos_file = "#{temp_dir}/mosdepth_frag.per-base.bed.gz"
      mop_file = "#{temp_dir}/mopdepth_frag.per-base.bed.gz"
      File.exists?(mos_file).should be_true
      File.exists?(mop_file).should be_true

      # Compare contents
      compare_files(mos_file, mop_file).should be_true
    end
  end

  describe "error handling" do
    it "handles missing chromosome gracefully" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["-c", "nonexistent"], "#{temp_dir}/mosdepth_missing", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["-c", "nonexistent"], "#{temp_dir}/mopdepth_missing", temp_dir, test_bam)

      # Both should fail with non-zero exit code
      mosdepth_status.success?.should be_false
      mopdepth_status.success?.should be_false
    end

    it "handles invalid arguments" do
      next unless mosdepth_available?(mosdepth_dir)
      mosdepth_status = run_mosdepth(["--invalid-option"], "#{temp_dir}/mosdepth_invalid", temp_dir, mosdepth_dir, test_bam)
      mopdepth_status = run_mopdepth(["--invalid-option"], "#{temp_dir}/mopdepth_invalid", temp_dir, test_bam)

      # Both should fail with non-zero exit code
      mosdepth_status.success?.should be_false
      mopdepth_status.success?.should be_false
    end
  end
end
