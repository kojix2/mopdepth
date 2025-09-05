require "./spec_helper"
require "file_utils"

describe "CSI index generation" do
  temp_dir = "/tmp/mopdepth_csi_test"
  test_bam = "./mosdepth/tests/ovl.bam"

  before_each do
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "creates .csi files by default" do
    next unless File.exists?(test_bam)
    status = TestBin.run(["-M"], "#{temp_dir}/csi_out", test_bam, temp_dir)
    status.success?.should be_true

    [
      "#{temp_dir}/csi_out.per-base.bed.gz.csi",
      "#{temp_dir}/csi_out.mosdepth.summary.txt",
    ].each do |path|
      # summary file must exist; .csi should exist for per-base
      if path.ends_with?(".csi")
        File.exists?(path).should be_true
      else
        File.exists?(path).should be_true
      end
    end
  end
end
