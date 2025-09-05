require "./spec_helper"
require "../src/depth/config"

describe Depth::Config do
  describe "initialization" do
    it "initializes with default values" do
      config = Depth::Config.new
      config.prefix.should eq("")
      config.path.should eq("")
      config.threads.should eq(0)
      config.chrom.should eq("")
      config.by.should eq("")
      config.no_per_base?.should be_false
      config.mapq.should eq(0)
      config.min_frag_len.should eq(-1)
      config.max_frag_len.should eq(-1)
      config.fast_mode?.should be_false
      config.fragment_mode?.should be_false
      config.use_median?.should be_false
      config.thresholds.should eq([] of Int32)
      config.quantize.should eq("")
    end
  end

  describe "#validate!" do
    it "raises error for empty path" do
      config = Depth::Config.new
      config.prefix = "test"

      expect_raises(ArgumentError, "BAM/CRAM path is required") do
        config.validate!
      end
    end

    it "raises error for empty prefix" do
      config = Depth::Config.new
      config.path = "test.bam"

      expect_raises(ArgumentError, "Output prefix is required") do
        config.validate!
      end
    end

    it "raises error for negative MAPQ" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.mapq = -1

      expect_raises(ArgumentError, "Invalid MAPQ threshold") do
        config.validate!
      end
    end

    it "raises error for negative threads" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.threads = -1

      expect_raises(ArgumentError, "Invalid thread count") do
        config.validate!
      end
    end

    it "raises error when min_frag_len > max_frag_len" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.min_frag_len = 100
      config.max_frag_len = 50

      expect_raises(ArgumentError, "min_frag_len cannot be greater than max_frag_len") do
        config.validate!
      end
    end

    it "raises error when both fast_mode and fragment_mode are true" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.fast_mode = true
      config.fragment_mode = true

      expect_raises(ArgumentError, "--fast-mode and --fragment-mode cannot be used together") do
        config.validate!
      end
    end

    it "validates successfully with valid configuration" do
      config = Depth::Config.new
      config.prefix = "test"
      config.path = "test.bam"
      config.mapq = 10
      config.threads = 4
      config.min_frag_len = 50
      config.max_frag_len = 500

      config.validate! # Should not raise
    end
  end

  describe "#has_regions?" do
    it "returns false for empty by" do
      config = Depth::Config.new
      config.has_regions?.should be_false
    end

    it "returns true for non-empty by" do
      config = Depth::Config.new
      config.by = "1000"
      config.has_regions?.should be_true
    end
  end

  describe "#window_size" do
    it "returns 0 for empty by" do
      config = Depth::Config.new
      config.window_size.should eq(0)
    end

    it "returns numeric value for numeric by" do
      config = Depth::Config.new
      config.by = "1000"
      config.window_size.should eq(1000)
    end

    it "returns 0 for non-numeric by" do
      config = Depth::Config.new
      config.by = "regions.bed"
      config.window_size.should eq(0)
    end
  end

  describe "#bed_path" do
    it "returns nil for empty by" do
      config = Depth::Config.new
      config.bed_path.should be_nil
    end

    it "returns nil for numeric by" do
      config = Depth::Config.new
      config.by = "1000"
      config.bed_path.should be_nil
    end

    it "returns path for non-numeric by" do
      config = Depth::Config.new
      config.by = "regions.bed"
      config.bed_path.should eq("regions.bed")
    end
  end

  describe "#has_quantize?" do
    it "returns false for empty quantize" do
      config = Depth::Config.new
      config.has_quantize?.should be_false
    end

    it "returns false for 'nil' quantize" do
      config = Depth::Config.new
      config.quantize = "nil"
      config.has_quantize?.should be_false
    end

    it "returns true for valid quantize string" do
      config = Depth::Config.new
      config.quantize = "0:1:4:"
      config.has_quantize?.should be_true
    end
  end

  describe "#quantize_args" do
    it "returns empty array for no quantize" do
      config = Depth::Config.new
      config.quantize_args.should eq([] of Int32)
    end

    it "returns parsed quantize args" do
      config = Depth::Config.new
      config.quantize = "0:1:4:"
      config.quantize_args.should eq([0, 1, 4, Int32::MAX])
    end

    it "handles simple quantize string" do
      config = Depth::Config.new
      config.quantize = ":10"
      config.quantize_args.should eq([0, 10])
    end
  end

  describe "#to_options" do
    it "converts configuration to options correctly" do
      config = Depth::Config.new
      config.mapq = 20
      config.min_frag_len = 100
      config.max_frag_len = 500
      config.fast_mode = true
      config.fragment_mode = false

      options = config.to_options
      options.mapq.should eq(20)
      options.min_frag_len.should eq(100)
      options.max_frag_len.should eq(500)
      options.fast_mode.should be_true
      options.fragment_mode.should be_false
    end

    it "handles negative max_frag_len" do
      config = Depth::Config.new
      config.max_frag_len = -1

      options = config.to_options
      options.max_frag_len.should eq(Int32::MAX)
    end
  end

  describe "integration with quantize functionality" do
    it "integrates with quantize module correctly" do
      config = Depth::Config.new
      config.quantize = "0:1:4:"

      config.has_quantize?.should be_true
      args = config.quantize_args
      args.should eq([0, 1, 4, Int32::MAX])

      # Test that it can create lookup table
      lookup = Depth::Stats::Quantize.make_lookup(args)
      lookup.size.should eq(3)
      lookup[2].should eq("4:inf")
    end
  end
end
