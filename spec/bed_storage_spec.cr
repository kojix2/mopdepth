require "./spec_helper"
require "../src/depth/io/bed_reader"

describe "BED interval storage" do
  it "keeps names, overlaps, and sorted coordinates with compact intervals" do
    path = File.tempname("mopdepth-bed-storage", ".bed")
    begin
      File.write(path, "chr1\t5\t9\tsecond\textra\tignored\nchr2\t0\t1\nchr1\t2\t7\tfirst\n")
      targets = [Depth::Core::Target.new("chr1", 10, 0), Depth::Core::Target.new("chr2", 10, 1)]
      regions = Depth::FileIO.read_bed(path, targets)
      regions["chr1"].map { |region| {region.start, region.stop, region.name} }.should eq([
        {2, 7, "first"}, {5, 9, "second"},
      ])
      regions["chr2"].first.name.should be_nil
    ensure
      File.delete(path) if File.exists?(path)
    end
  end
end
