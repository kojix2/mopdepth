require "./spec_helper"
require "../src/depth/stats/depth_stat"

describe Depth::Stats::DepthStat do
  describe "initialization" do
    it "initializes with correct default values" do
      stat = Depth::Stats::DepthStat.new
      stat.n_bases.should eq(0)
      stat.sum_depth.should eq(0_i64)
      stat.min_depth.should eq(Int32::MAX)
      stat.max_depth.should eq(0)
    end
  end

  describe "#clear" do
    it "resets all values to defaults" do
      stat = Depth::Stats::DepthStat.new
      stat.n_bases = 100
      stat.sum_depth = 500_i64
      stat.min_depth = 1
      stat.max_depth = 10

      stat.clear

      stat.n_bases.should eq(0)
      stat.sum_depth.should eq(0_i64)
      stat.min_depth.should eq(Int32::MAX)
      stat.max_depth.should eq(0)
    end

    it "sets min_depth to maximum value after clear (like original mosdepth)" do
      stat = Depth::Stats::DepthStat.new
      stat.min_depth = 5

      stat.clear

      stat.min_depth.should be > 0
      stat.min_depth.should eq(Int32::MAX)
    end
  end

  describe ".from_array" do
    it "creates correct statistics from empty array" do
      data = [] of Int32
      stat = Depth::Stats::DepthStat.from_array(data)

      stat.n_bases.should eq(0)
      stat.sum_depth.should eq(0_i64)
      stat.min_depth.should eq(Int32::MAX)
      stat.max_depth.should eq(0)
    end

    it "creates correct statistics from array with single element" do
      data = [5]
      stat = Depth::Stats::DepthStat.from_array(data)

      stat.n_bases.should eq(1)
      stat.sum_depth.should eq(5_i64)
      stat.min_depth.should eq(5)
      stat.max_depth.should eq(5)
    end

    it "creates correct statistics from array with multiple elements" do
      data = [1, 5, 3, 8, 2]
      stat = Depth::Stats::DepthStat.from_array(data)

      stat.n_bases.should eq(5)
      stat.sum_depth.should eq(19_i64)
      stat.min_depth.should eq(1)
      stat.max_depth.should eq(8)
    end

    it "handles array with zeros correctly" do
      data = [0, 5, 0, 3, 0]
      stat = Depth::Stats::DepthStat.from_array(data)

      stat.n_bases.should eq(5)
      stat.sum_depth.should eq(8_i64)
      stat.min_depth.should eq(0)
      stat.max_depth.should eq(5)
    end

    it "works with start and end indices" do
      data = [1, 2, 3, 4, 5]
      stat = Depth::Stats::DepthStat.from_array(data, 1, 3)

      stat.n_bases.should eq(3)       # indices 1, 2, 3
      stat.sum_depth.should eq(9_i64) # 2 + 3 + 4
      stat.min_depth.should eq(2)
      stat.max_depth.should eq(4)
    end

    it "handles negative end index (uses array size - 1)" do
      data = [1, 2, 3, 4, 5]
      stat = Depth::Stats::DepthStat.from_array(data, 0, -1)

      stat.n_bases.should eq(5)
      stat.sum_depth.should eq(15_i64)
      stat.min_depth.should eq(1)
      stat.max_depth.should eq(5)
    end
  end

  describe ".from_slice" do
    it "creates correct statistics from slice" do
      data = [1, 5, 3, 8, 2]
      slice = Slice.new(data.to_unsafe, data.size)
      stat = Depth::Stats::DepthStat.from_slice(slice)

      stat.n_bases.should eq(5)
      stat.sum_depth.should eq(19_i64)
      stat.min_depth.should eq(1)
      stat.max_depth.should eq(8)
    end

    it "creates correct statistics from empty slice" do
      data = [] of Int32
      slice = Slice.new(data.to_unsafe, data.size)
      stat = Depth::Stats::DepthStat.from_slice(slice)

      stat.n_bases.should eq(0)
      stat.sum_depth.should eq(0_i64)
      stat.min_depth.should eq(Int32::MAX)
      stat.max_depth.should eq(0)
    end
  end

  describe "#+" do
    it "combines two DepthStats correctly" do
      stat1 = Depth::Stats::DepthStat.new
      stat1.n_bases = 3
      stat1.sum_depth = 10_i64
      stat1.min_depth = 1
      stat1.max_depth = 5

      stat2 = Depth::Stats::DepthStat.new
      stat2.n_bases = 2
      stat2.sum_depth = 8_i64
      stat2.min_depth = 2
      stat2.max_depth = 6

      result = stat1 + stat2

      result.n_bases.should eq(5)
      result.sum_depth.should eq(18_i64)
      result.min_depth.should eq(1) # min of 1 and 2
      result.max_depth.should eq(6) # max of 5 and 6
    end

    it "handles combining with empty stat" do
      stat1 = Depth::Stats::DepthStat.from_array([1, 2, 3])
      stat2 = Depth::Stats::DepthStat.new # empty

      result = stat1 + stat2

      result.n_bases.should eq(3)
      result.sum_depth.should eq(6_i64)
      result.min_depth.should eq(1)
      result.max_depth.should eq(3)
    end

    it "handles edge case with Int32::MAX min_depth" do
      stat1 = Depth::Stats::DepthStat.new
      stat1.n_bases = 2
      stat1.sum_depth = 5_i64
      stat1.min_depth = 2
      stat1.max_depth = 3

      stat2 = Depth::Stats::DepthStat.new # min_depth is Int32::MAX

      result = stat1 + stat2

      result.min_depth.should eq(2) # should take the smaller value
    end
  end

  describe "integration with original mosdepth behavior" do
    it "behaves like original mosdepth depthstat min test" do
      # Test equivalent to: var d = newSeq[int32](); var t = newDepthStat(d); check t.min_depth > 0
      data = [] of Int32
      stat = Depth::Stats::DepthStat.from_array(data)
      stat.min_depth.should be > 0
      stat.min_depth.should eq(Int32::MAX)
    end

    it "behaves like original mosdepth clear test" do
      # Test equivalent to: var dd: depth_stat; check dd.min_depth == 0; dd.clear(); check dd.min_depth > 0
      stat = Depth::Stats::DepthStat.new
      # In Crystal, initial min_depth is Int32::MAX, not 0 like in Nim
      # But after clear(), it should still be > 0
      stat.clear
      stat.min_depth.should be > 0
      stat.min_depth.should eq(Int32::MAX)
    end
  end
end
