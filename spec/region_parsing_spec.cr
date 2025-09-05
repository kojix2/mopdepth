require "./spec_helper"
require "../src/depth/io/bed_reader"

describe "Region Parsing" do
  describe "parse_region_str" do
    it "parses chromosome only" do
      region = Depth::FileIO.parse_region_str("Super-Scaffold_52")
      region.should_not be_nil
      if region
        region.chrom.should eq("Super-Scaffold_52")
        region.start.should eq(0)
        region.stop.should eq(0)
      end
    end

    it "parses chromosome with range" do
      region = Depth::FileIO.parse_region_str("Super-Scaffold_52:2-1000")
      if region
        region.chrom.should eq("Super-Scaffold_52")
        region.start.should eq(1) # 1-based to 0-based conversion
        region.stop.should eq(1000)
      end
    end

    it "parses chromosome with single position" do
      region = Depth::FileIO.parse_region_str("chr1:100")
      region.should_not be_nil
      if region
        region.chrom.should eq("chr1")
        region.start.should eq(99) # 1-based to 0-based conversion
        region.stop.should eq(100)
      end
    end

    it "returns nil for empty string" do
      region = Depth::FileIO.parse_region_str("")
      region.should be_nil
    end

    it "returns nil for nil string" do
      region = Depth::FileIO.parse_region_str("nil")
      region.should be_nil
    end

    it "handles complex chromosome names" do
      region = Depth::FileIO.parse_region_str("chr22:20000000-23000000")
      region.should_not be_nil
      if region
        region.chrom.should eq("chr22")
        region.start.should eq(19999999) # 1-based to 0-based conversion
        region.stop.should eq(23000000)
      end
    end

    it "handles mitochondrial chromosome" do
      region = Depth::FileIO.parse_region_str("MT:1-16569")
      region.should_not be_nil
      if region
        region.chrom.should eq("MT")
        region.start.should eq(0) # 1-based to 0-based conversion
        region.stop.should eq(16569)
      end
    end
  end
end
