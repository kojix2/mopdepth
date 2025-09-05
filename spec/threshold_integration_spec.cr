require "./spec_helper"
require "../src/depth/config"
require "../src/depth/io/output_manager"
require "../src/depth/stats/threshold"

describe "Threshold Integration" do
  describe "Config integration" do
    it "correctly identifies when thresholds are enabled" do
      config = Depth::Config.new
      config.thresholds_str = "1,2,3"
      config.has_thresholds?.should be_true
    end

    it "correctly identifies when thresholds are disabled" do
      config = Depth::Config.new
      config.thresholds_str = ""
      config.has_thresholds?.should be_false
    end

    it "returns correct threshold values" do
      config = Depth::Config.new
      config.thresholds_str = "1,5,10"
      config.threshold_values.should eq([1, 5, 10])
    end

    it "validates threshold requires regions" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.thresholds_str = "1,2,3"
      config.by = "" # No regions specified

      expect_raises(ArgumentError, "--thresholds can only be used when --by is specified") do
        config.validate!
      end
    end

    it "validates successfully with thresholds and regions" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.thresholds_str = "1,2,3"
      config.by = "100" # Window size specified

      config.validate! # Should not raise
    end
  end

  describe "OutputManager with thresholds" do
    it "creates threshold output file when has_thresholds is true" do
      prefix = "test_threshold"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.by = "100" # regions required for thresholds
      config.thresholds_str = "1,5,10"

      output = Depth::FileIO::OutputManager.new(config)

      begin
        output.f_thresholds.should_not be_nil
        File.exists?("#{prefix}.thresholds.bed.gz").should be_true
      ensure
        output.close_all
        File.delete("#{prefix}.thresholds.bed.gz") if File.exists?("#{prefix}.thresholds.bed.gz")
        # remove both labels (mopdepth and mosdepth)
        [
          "#{prefix}.mopdepth.summary.txt",
          "#{prefix}.mopdepth.global.dist.txt",
          "#{prefix}.mopdepth.region.dist.txt",
          "#{prefix}.mosdepth.summary.txt",
          "#{prefix}.mosdepth.global.dist.txt",
          "#{prefix}.mosdepth.region.dist.txt",
        ].each { |f| File.delete(f) if File.exists?(f) }
        File.delete("#{prefix}.regions.bed.gz") if File.exists?("#{prefix}.regions.bed.gz")
        File.delete("#{prefix}.per-base.bed.gz") if File.exists?("#{prefix}.per-base.bed.gz")
      end
    end

    it "does not create threshold output file when has_thresholds is false" do
      prefix = "test_no_threshold"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.by = "100" # regions specified but no thresholds

      output = Depth::FileIO::OutputManager.new(config)

      begin
        output.f_thresholds.should be_nil
        File.exists?("#{prefix}.thresholds.bed.gz").should be_false
      ensure
        output.close_all
        [
          "#{prefix}.mopdepth.summary.txt",
          "#{prefix}.mopdepth.global.dist.txt",
          "#{prefix}.mopdepth.region.dist.txt",
          "#{prefix}.mosdepth.summary.txt",
          "#{prefix}.mosdepth.global.dist.txt",
          "#{prefix}.mosdepth.region.dist.txt",
        ].each { |f| File.delete(f) if File.exists?(f) }
        File.delete("#{prefix}.regions.bed.gz") if File.exists?("#{prefix}.regions.bed.gz")
        File.delete("#{prefix}.per-base.bed.gz") if File.exists?("#{prefix}.per-base.bed.gz")
      end
    end

    it "writes threshold header correctly" do
      prefix = "test_header"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.by = "100"
      config.thresholds_str = "1,5,10"

      output = Depth::FileIO::OutputManager.new(config)

      begin
        thresholds = [1, 5, 10]
        output.write_thresholds_header(thresholds)
        output.close_all

        content = TestIO.read_text("#{prefix}.thresholds.bed.gz")
        content.should eq("#chrom\tstart\tend\tregion\t1X\t5X\t10X\n")
      ensure
        File.delete("#{prefix}.thresholds.bed.gz") if File.exists?("#{prefix}.thresholds.bed.gz")
        [
          "#{prefix}.mopdepth.summary.txt",
          "#{prefix}.mopdepth.global.dist.txt",
          "#{prefix}.mopdepth.region.dist.txt",
          "#{prefix}.mosdepth.summary.txt",
          "#{prefix}.mosdepth.global.dist.txt",
          "#{prefix}.mosdepth.region.dist.txt",
        ].each { |f| File.delete(f) if File.exists?(f) }
        File.delete("#{prefix}.regions.bed.gz") if File.exists?("#{prefix}.regions.bed.gz")
        File.delete("#{prefix}.per-base.bed.gz") if File.exists?("#{prefix}.per-base.bed.gz")
      end
    end

    it "writes threshold counts correctly" do
      prefix = "test_counts"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.by = "100"
      config.thresholds_str = "1,5,10"

      output = Depth::FileIO::OutputManager.new(config)

      begin
        thresholds = [1, 5, 10]
        output.write_thresholds_header(thresholds)

        counts = [100, 80, 20]
        output.write_threshold_counts("chr1", 0, 100, "region1", counts)
        output.write_threshold_counts("chr1", 100, 200, nil, [50, 30, 10])
        output.close_all

        content = TestIO.read_text("#{prefix}.thresholds.bed.gz")
        lines = content.split('\n')
        lines[0].should eq("#chrom\tstart\tend\tregion\t1X\t5X\t10X")
        lines[1].should eq("chr1\t0\t100\tregion1\t100\t80\t20")
        lines[2].should eq("chr1\t100\t200\tunknown\t50\t30\t10")
      ensure
        File.delete("#{prefix}.thresholds.bed.gz") if File.exists?("#{prefix}.thresholds.bed.gz")
        File.delete("#{prefix}.mopdepth.summary.txt") if File.exists?("#{prefix}.mopdepth.summary.txt")
        File.delete("#{prefix}.mopdepth.global.dist.txt") if File.exists?("#{prefix}.mopdepth.global.dist.txt")
        File.delete("#{prefix}.mopdepth.region.dist.txt") if File.exists?("#{prefix}.mopdepth.region.dist.txt")
        File.delete("#{prefix}.regions.bed.gz") if File.exists?("#{prefix}.regions.bed.gz")
        File.delete("#{prefix}.per-base.bed.gz") if File.exists?("#{prefix}.per-base.bed.gz")
      end
    end
  end

  describe "Compatibility with original mosdepth" do
    it "produces same threshold args as original mosdepth test cases" do
      # Test case from mosdepth: threshold_args("1,2,3")
      ts = Depth::Stats.threshold_args("1,2,3")
      ts[0].should eq(1)
      ts[1].should eq(2)
      ts[2].should eq(3)
      ts.size.should eq(3)
    end

    it "handles various threshold configurations like original mosdepth" do
      # Empty case
      ts = Depth::Stats.threshold_args("")
      ts.should eq([] of Int32)

      # Single value
      ts = Depth::Stats.threshold_args("5")
      ts.should eq([5])

      # Multiple values with sorting
      ts = Depth::Stats.threshold_args("10,1,5")
      ts.should eq([1, 5, 10])
    end
  end
end
