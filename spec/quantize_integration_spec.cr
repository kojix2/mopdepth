require "./spec_helper"
require "../src/depth/config"
require "../src/depth/io/output_manager"
require "../src/depth/stats/quantize"

# Helper module for test file cleanup
module TestCleanup
  def self.cleanup_test_files(prefix : String)
    files_to_clean = [
      "#{prefix}.quantized.bed.gz",
      "#{prefix}.per-base.bed.gz",
      # summary/dist both styles
      "#{prefix}.mopdepth.summary.txt",
      "#{prefix}.mopdepth.global.dist.txt",
      "#{prefix}.mopdepth.region.dist.txt",
      "#{prefix}.mosdepth.summary.txt",
      "#{prefix}.mosdepth.global.dist.txt",
      "#{prefix}.mosdepth.region.dist.txt",
      "#{prefix}.regions.bed.gz",
      "#{prefix}.thresholds.bed.gz",
    ]

    files_to_clean.each do |file|
      File.delete(file) if File.exists?(file)
    end
  end
end

describe "Quantize Integration" do
  after_each do
    # Clean up any test files that might have been created
    TestCleanup.cleanup_test_files("test_quantize")
    TestCleanup.cleanup_test_files("test_no_quantize")
  end
  describe "OutputManager with quantize" do
    it "creates quantized output file when has_quantize is true" do
      prefix = "test_quantize"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.quantize = "0:1:4:"

      output = Depth::FileIO::OutputManager.new(config)

      # Check that quantized file was created
      output.f_quantized.should_not be_nil

      # Write some test data
      output.write_quantized_interval("chr1", 0, 100, "0:1")
      output.write_quantized_interval("chr1", 100, 200, "1:4")

      output.close_all

      # Verify file exists and has content
      File.exists?("#{prefix}.quantized.bed.gz").should be_true
      content = TestIO.read_text("#{prefix}.quantized.bed.gz")
      content.should contain("chr1\t0\t100\t0:1")
      content.should contain("chr1\t100\t200\t1:4")
    end

    it "does not create quantized output file when has_quantize is false" do
      prefix = "test_no_quantize"
      config = Depth::Config.new
      config.prefix = prefix
      config.path = "test.bam"
      config.quantize = ""

      output = Depth::FileIO::OutputManager.new(config)

      # Check that quantized file was not created
      output.f_quantized.should be_nil

      output.close_all

      # Verify file does not exist
      File.exists?("#{prefix}.quantized.bed.gz").should be_false
    end
  end

  describe "Config integration" do
    it "correctly identifies when quantize is enabled" do
      config = Depth::Config.new
      config.quantize = "0:1:4:"

      config.has_quantize?.should be_true
      config.quantize_args.should eq([0, 1, 4, Int32::MAX])
    end

    it "correctly identifies when quantize is disabled" do
      config = Depth::Config.new
      config.quantize = ""

      config.has_quantize?.should be_false
      config.quantize_args.should eq([] of Int32)
    end
  end

  describe "End-to-end quantize workflow" do
    it "processes quantize workflow correctly" do
      # Test the complete workflow from configuration to output
      config = Depth::Config.new
      config.quantize = "0:1:4:"

      # Get quantize args
      quants = config.quantize_args
      quants.should eq([0, 1, 4, Int32::MAX])

      # Create lookup table
      lookup = Depth::Stats::Quantize.make_lookup(quants)
      lookup.should eq(["0:1", "1:4", "4:inf"])

      # Test quantized generation with sample coverage data
      coverage = [0, 0, 1, 1, 1, 4, 4, 0, 0]
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      # Should generate segments based on quantization changes
      segments.size.should be > 0
      segments.each do |start, stop, label|
        start.should be >= 0
        stop.should be > start
        label.should match(/\d+:(inf|\d+)/)
      end
    end

    it "handles edge case with no data correctly" do
      config = Depth::Config.new
      config.quantize = "0:1:4:"

      quants = config.quantize_args

      # Test with empty coverage (simulating no data case)
      coverage = [] of Int32
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      # Should handle empty coverage gracefully
      segments.should be_empty
    end

    it "respects environment variables for custom labels" do
      # Set environment variable
      ENV["MOSDEPTH_Q0"] = "LOW"
      ENV["MOSDEPTH_Q1"] = "MEDIUM"
      ENV["MOSDEPTH_Q2"] = "HIGH"

      begin
        quants = [0, 1, 4, Int32::MAX]
        lookup = Depth::Stats::Quantize.make_lookup(quants)

        lookup[0].should eq("LOW")
        lookup[1].should eq("MEDIUM")
        lookup[2].should eq("HIGH")
      ensure
        # Clean up environment variables
        ENV.delete("MOSDEPTH_Q0")
        ENV.delete("MOSDEPTH_Q1")
        ENV.delete("MOSDEPTH_Q2")
      end
    end
  end

  describe "Compatibility with original mosdepth" do
    it "produces same quantize args as original mosdepth test cases" do
      # Test cases from original mosdepth functional-tests.sh

      # Test case: -q 0:1:1000
      config = Depth::Config.new
      config.quantize = "0:1:1000"
      args = config.quantize_args
      args.should eq([0, 1, 1000])

      lookup = Depth::Stats::Quantize.make_lookup(args)
      lookup.should eq(["0:1", "1:1000"])
    end

    it "handles single quantize value like original mosdepth" do
      # Test case: -q 60 (single threshold)
      config = Depth::Config.new
      config.quantize = "60"
      args = config.quantize_args
      args.should eq([0, 60, Int32::MAX])

      lookup = Depth::Stats::Quantize.make_lookup(args)
      lookup.should eq(["0:60", "60:inf"])
    end
  end
end
