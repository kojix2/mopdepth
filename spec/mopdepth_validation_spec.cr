require "./spec_helper"
require "file_utils"
require "process"

# Minimal smoke tests: only assert that expected output files (gz) are created
def run_mopdepth_min(args : Array(String), prefix : String, bam : String, temp_dir : String) : Process::Status
  TestBin.run(args, prefix, bam, temp_dir)
end

describe "mopdepth minimal output smoke tests (gz only, no -M)" do
  temp_dir = "/tmp/mopdepth_validation"
  mosdepth_dir = "./mosdepth"
  test_bam = "#{mosdepth_dir}/tests/ovl.bam"

  before_each do
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "writes per-base .bed.gz by default" do
    status = run_mopdepth_min([] of String, "#{temp_dir}/test", test_bam, temp_dir)
    status.success?.should be_true
    File.exists?("#{temp_dir}/test.per-base.bed.gz").should be_true
  end

  it "writes summary and global dist (mopdepth.* label without -M)" do
    status = run_mopdepth_min([] of String, "#{temp_dir}/sumdist", test_bam, temp_dir)
    status.success?.should be_true
    File.exists?("#{temp_dir}/sumdist.mopdepth.summary.txt").should be_true
    File.exists?("#{temp_dir}/sumdist.mopdepth.global.dist.txt").should be_true
  end

  it "writes regions .bed.gz with -b windows" do
    status = run_mopdepth_min(["-b", "100"], "#{temp_dir}/win", test_bam, temp_dir)
    status.success?.should be_true
    File.exists?("#{temp_dir}/win.regions.bed.gz").should be_true
  end

  it "writes quantized .bed.gz with -q" do
    status = run_mopdepth_min(["-q", "0:1:1000"], "#{temp_dir}/qz", test_bam, temp_dir)
    status.success?.should be_true
    File.exists?("#{temp_dir}/qz.quantized.bed.gz").should be_true
  end

  it "writes thresholds .bed.gz with -T and -b" do
    status = run_mopdepth_min(["-T", "0,1,2", "-b", "100"], "#{temp_dir}/th", test_bam, temp_dir)
    status.success?.should be_true
    File.exists?("#{temp_dir}/th.thresholds.bed.gz").should be_true
  end
end
