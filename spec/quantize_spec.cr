require "./spec_helper"
require "../src/depth/stats/quantize"

describe Depth::Stats::Quantize do
  describe ".get_quantize_args" do
    it "parses simple quantize string" do
      result = Depth::Stats::Quantize.get_quantize_args(":1")
      result.should eq([0, 1])
    end

    it "parses complex quantize string" do
      result = Depth::Stats::Quantize.get_quantize_args("0:1:4:")
      result.should eq([0, 1, 4, Int32::MAX])
    end

    it "handles single number by wrapping with colons" do
      result = Depth::Stats::Quantize.get_quantize_args("5")
      result.should eq([0, 5, Int32::MAX])
    end

    it "handles string starting with colon" do
      result = Depth::Stats::Quantize.get_quantize_args(":10:")
      result.should eq([0, 10, Int32::MAX])
    end

    it "handles string ending with colon" do
      result = Depth::Stats::Quantize.get_quantize_args("5:10:")
      result.should eq([5, 10, Int32::MAX])
    end

    it "sorts the values" do
      result = Depth::Stats::Quantize.get_quantize_args("10:1:5")
      result.should eq([1, 5, 10])
    end

    it "returns empty array for nil or empty string" do
      Depth::Stats::Quantize.get_quantize_args("nil").should eq([] of Int32)
      Depth::Stats::Quantize.get_quantize_args("").should eq([] of Int32)
    end

    it "handles duplicate values" do
      result = Depth::Stats::Quantize.get_quantize_args("1:1:5")
      result.should eq([1, 1, 5])
    end
  end

  describe ".linear_search" do
    it "finds correct bin for values within range" do
      bins = [10, 22, 44, 99]

      Depth::Stats::Quantize.linear_search(10, bins).should eq(0)
      Depth::Stats::Quantize.linear_search(22, bins).should eq(1)
      Depth::Stats::Quantize.linear_search(44, bins).should eq(2)
      Depth::Stats::Quantize.linear_search(99, bins).should eq(3)
    end

    it "returns -1 for values outside range" do
      bins = [10, 22, 44, 99]

      Depth::Stats::Quantize.linear_search(8, bins).should eq(-1)
      Depth::Stats::Quantize.linear_search(800, bins).should eq(-1)
    end

    it "handles edge cases with special bins" do
      bins = [0, 1, Int32::MAX]

      Depth::Stats::Quantize.linear_search(0, bins).should eq(0)
      Depth::Stats::Quantize.linear_search(-1, bins).should eq(-1)
      Depth::Stats::Quantize.linear_search(99999, bins).should eq(1)
    end

    it "returns -1 for empty array" do
      Depth::Stats::Quantize.linear_search(5, [] of Int32).should eq(-1)
    end

    it "finds correct bin for values between boundaries" do
      bins = [10, 22, 44, 99]

      Depth::Stats::Quantize.linear_search(15, bins).should eq(0) # between 10 and 22
      Depth::Stats::Quantize.linear_search(30, bins).should eq(1) # between 22 and 44
      Depth::Stats::Quantize.linear_search(50, bins).should eq(2) # between 44 and 99
    end
  end

  describe ".make_lookup" do
    it "creates correct lookup table" do
      bins = [10, 22, 44, 99]
      lookup = Depth::Stats::Quantize.make_lookup(bins)

      lookup.should eq(["10:22", "22:44", "44:99"])
    end

    it "handles bins with infinity" do
      bins = [0, 10]
      lookup = Depth::Stats::Quantize.make_lookup(bins)

      lookup.should eq(["0:10"])
    end

    it "handles bins ending with Int32::MAX" do
      bins = Depth::Stats::Quantize.get_quantize_args("0:1:4:")
      lookup = Depth::Stats::Quantize.make_lookup(bins)

      lookup.size.should eq(3)
      lookup[2].should eq("4:inf")
    end

    it "returns empty array for single element" do
      bins = [10]
      lookup = Depth::Stats::Quantize.make_lookup(bins)

      lookup.should eq([] of String)
    end
  end

  describe ".gen_quantized" do
    it "generates quantized segments correctly" do
      quants = [0, 1, 4]
      coverage = [0, 0, 1, 1, 1, 4, 4, 0, 0]
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      # Should generate segments based on changes in quantization bins
      segments.size.should be > 0
      segments.each do |start, stop, label|
        start.should be >= 0
        stop.should be > start
        label.should_not be_empty
      end
    end

    it "handles empty coverage array" do
      quants = [0, 1, 4]
      coverage = [] of Int32
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      segments.should be_empty
    end

    it "handles empty quantization array" do
      quants = [] of Int32
      coverage = [1, 2, 3]
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      segments.should be_empty
    end

    it "generates correct labels" do
      quants = [0, 1]
      coverage = [0, 1, 0]
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      segments.each do |_start, _stop, label|
        label.should match(/\d+:(inf|\d+)/)
      end
    end
  end

  describe "integration with original mosdepth behavior" do
    it "behaves like original mosdepth quantize-args test" do
      # Test equivalent to original mosdepth test cases
      rs = Depth::Stats::Quantize.get_quantize_args(":1")
      rs[0].should eq(0)
      rs[1].should eq(1)

      rs = Depth::Stats::Quantize.get_quantize_args("0:1:4:")
      rs[0].should eq(0)
      rs[1].should eq(1)
      rs[2].should eq(4)
      rs[3].should eq(Int32::MAX)
    end

    it "behaves like original mosdepth linear-search test" do
      bins = [10, 22, 44, 99]

      bins.each_with_index do |v, i|
        idx = Depth::Stats::Quantize.linear_search(v, bins)
        idx.should eq(i)
      end

      Depth::Stats::Quantize.linear_search(8, bins).should eq(-1)
      Depth::Stats::Quantize.linear_search(800, bins).should eq(-1)

      bins = [0, 1, Int32::MAX]
      Depth::Stats::Quantize.linear_search(0, bins).should eq(0)
      Depth::Stats::Quantize.linear_search(-1, bins).should eq(-1)
      Depth::Stats::Quantize.linear_search(99999, bins).should eq(1)
    end

    it "behaves like original mosdepth lookup test" do
      bins = [10, 22, 44, 99]
      lookup = Depth::Stats::Quantize.make_lookup(bins)
      lookup[0].should eq("10:22")
      lookup[1].should eq("22:44")
      lookup[2].should eq("44:99")
      lookup.size.should eq(3)

      bins = [0, 10]
      lookup = Depth::Stats::Quantize.make_lookup(bins)
      lookup[0].should eq("0:10")
      lookup.size.should eq(1)

      bins = Depth::Stats::Quantize.get_quantize_args("0:1:4:")
      lookup = Depth::Stats::Quantize.make_lookup(bins)
      lookup.size.should eq(3)
      lookup[2].should eq("4:inf")
    end
  end
end
