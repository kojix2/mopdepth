require "./spec_helper"
require "../src/mopdepth"

private def create_library_bam(path : String)
  header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:empty\tLN:5\n")
  HTS::Bam.open(path, "wb") do |bam|
    bam.write_header(header)
    bam.write(HTS::Bam::Record.new(
      header, "read", 0_u16, "chr1", 2_i64, 60_u8, "3M", "AAA",
      [30_u8, 30_u8, 30_u8], -1, 0_i64, 0_i64
    ))
  end
  HTS::Bam.build_index(path, verbose: false)
end

describe Mopdepth::Reader do
  temp_dir = "/tmp/mopdepth_library"

  before_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "exposes targets and owns the indexed input handle" do
    path = "#{temp_dir}/input.bam"
    create_library_bam(path)

    reader = Mopdepth::Reader.open(path)
    reader.targets.map { |target| {target.name, target.length} }.should eq([
      {"chr1", 10},
      {"empty", 5},
    ])
    reader.file_format.should eq("bam")
    reader.close
    reader.closed?.should be_true
    expect_raises(Mopdepth::ClosedError) { reader.coverage("chr1") }
  end

  it "returns dense coverage and reuses caller-owned storage" do
    path = "#{temp_dir}/coverage.bam"
    create_library_bam(path)
    progress = [] of Mopdepth::Progress
    control = Mopdepth::Control.new { |event| progress << event }

    Mopdepth::Reader.open(path) do |reader|
      buffer = [99_i32, 100_i32]
      result = reader.coverage_into(Mopdepth::Region.new("chr1", 1, 6), buffer, control: control)
      result.values.same?(buffer).should be_true
      result.values.should eq([0, 1, 1, 1, 0])
      result.observed?.should be_true
      result.records_examined.should eq(1)

      filtered = reader.coverage("chr1", options: Mopdepth::Options.new(minimum_mapping_quality: 61))
      filtered.values.should eq(Array(Int32).new(10, 0_i32))
      filtered.observed?.should be_true
    end

    progress.map(&.stage).should contain(Mopdepth::Stage::Preparing)
    progress.map(&.stage).should contain(Mopdepth::Stage::Reading)
    progress.map(&.stage).should contain(Mopdepth::Stage::Finalizing)
  end

  it "iterates references with independent result storage" do
    path = "#{temp_dir}/each.bam"
    create_library_bam(path)
    results = [] of Mopdepth::Coverage

    Mopdepth::Reader.open(path) do |reader|
      reader.each_coverage { |coverage| results << coverage }
    end

    results.map(&.region.chromosome).should eq(["chr1", "empty"])
    results[0].values.should eq([0, 0, 1, 1, 1, 0, 0, 0, 0, 0])
    results[1].values.should eq([0, 0, 0, 0, 0])
    results[1].status.no_data?.should be_true
    results[0].values.same?(results[1].values).should be_false
  end

  it "cancels cooperatively before querying records" do
    path = "#{temp_dir}/cancel.bam"
    create_library_bam(path)
    control = Mopdepth::Control.new
    control.cancel

    Mopdepth::Reader.open(path) do |reader|
      expect_raises(Mopdepth::CancelledError) do
        reader.coverage("chr1", control: control)
      end
    end
  end

  it "validates public options and regions" do
    expect_raises(ArgumentError) { Mopdepth::Options.new(minimum_mapping_quality: 256).validate! }
    expect_raises(ArgumentError) { Mopdepth::Options.new(minimum_fragment_length: 10, maximum_fragment_length: 9).validate! }
    expect_raises(ArgumentError) { Mopdepth::Region.new("chr1", -1, 2) }
  end
end
