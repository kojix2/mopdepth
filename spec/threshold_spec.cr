require "./spec_helper"
require "../src/depth/stats/threshold"

describe "Threshold" do
  describe "threshold_args" do
    it "parses comma-separated threshold values" do
      ts = Depth::Stats.threshold_args("1,2,3")
      ts.should eq([1, 2, 3])
      ts.size.should eq(3)
    end

    it "sorts threshold values" do
      ts = Depth::Stats.threshold_args("3,1,2")
      ts.should eq([1, 2, 3])
    end

    it "handles single threshold value" do
      ts = Depth::Stats.threshold_args("5")
      ts.should eq([5])
      ts.size.should eq(1)
    end

    it "returns empty array for empty string" do
      ts = Depth::Stats.threshold_args("")
      ts.should eq([] of Int32)
      ts.size.should eq(0)
    end

    it "returns empty array for nil string" do
      ts = Depth::Stats.threshold_args("nil")
      ts.should eq([] of Int32)
      ts.size.should eq(0)
    end

    it "handles duplicate values by keeping them" do
      ts = Depth::Stats.threshold_args("1,2,2,3")
      ts.should eq([1, 2, 2, 3])
    end

    it "handles larger threshold values" do
      ts = Depth::Stats.threshold_args("10,50,100")
      ts.should eq([10, 50, 100])
    end

    it "raises error for invalid threshold string" do
      expect_raises(ArgumentError, "Invalid threshold string: 'invalid'") do
        Depth::Stats.threshold_args("invalid")
      end
    end

    it "raises error for mixed valid and invalid values" do
      expect_raises(ArgumentError, "Invalid threshold string: '1,invalid,3'") do
        Depth::Stats.threshold_args("1,invalid,3")
      end
    end
  end
end
