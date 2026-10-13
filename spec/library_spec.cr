require "./spec_helper"
require "../src/mopdepth"

private def create_library_bam(path : String, records : Int32 = 1)
  header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:empty\tLN:5\n")
  HTS::Bam.open(path, "wb") do |bam|
    bam.write_header(header)
    record = HTS::Bam::Record.new(
      header, "read", 0_u16, "chr1", 2_i64, 60_u8, "3M", "AAA",
      [30_u8, 30_u8, 30_u8], -1, 0_i64, 0_i64
    )
    records.times { bam.write(record) }
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

  it "reuses a reader after cancellation inside an indexed query" do
    path = "#{temp_dir}/many.bam"
    create_library_bam(path, 20_000)
    controls = [] of Mopdepth::Control
    control = Mopdepth::Control.new do |event|
      controls.first.cancel if event.records_examined >= 8192
    end
    controls << control
    Mopdepth::Reader.open(path) do |reader|
      expect_raises(Mopdepth::CancelledError) { reader.coverage("chr1", control: control) }
      reader.records_examined.should eq(8192)
      reader.coverage("chr1").values.should eq([0, 0, 20_000, 20_000, 20_000, 0, 0, 0, 0, 0])
    end
  end

  it "resizes reused buffers and handles empty or invalid regions" do
    path = "#{temp_dir}/regions.bam"
    create_library_bam(path)
    Mopdepth::Reader.open(path) do |reader|
      values = reader.coverage("chr1").values
      reader.coverage_into(Mopdepth::Region.new("chr1", 3, 4), values).values.should eq([1])
      reader.coverage_into(reader.target("empty").region, values).values.should eq([0, 0, 0, 0, 0])
      reader.coverage_into(reader.target("chr1").region, values).values.should eq([0, 0, 1, 1, 1, 0, 0, 0, 0, 0])
      before = reader.records_examined
      reader.coverage_into(Mopdepth::Region.new("chr1", 10, 10), values).values.should be_empty
      reader.records_examined.should eq(before)
      expect_raises(Mopdepth::UnknownReferenceError) { reader.coverage("missing") }
      expect_raises(ArgumentError) { reader.coverage(Mopdepth::Region.new("chr1", 0, 11)) }
    end
  end

  it "counts a spanning fragment when neither mate overlaps the requested region" do
    path = "#{temp_dir}/fragments.bam"
    header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr\tLN:30\n")
    HTS::Bam.open(path, "wb") do |bam|
      bam.write_header(header)
      [{99_u16, 2_i64, 20_i64, 21_i64}, {147_u16, 20_i64, 2_i64, -21_i64}].each do |flag, position, mate, length|
        bam.write(HTS::Bam::Record.new(header, "pair", flag, "chr", position, 60_u8, "3M", "AAA",
          [30_u8] * 3, 0, mate, length))
      end
    end
    HTS::Bam.build_index(path, verbose: false)
    Mopdepth::Reader.open(path) do |reader|
      region = Mopdepth::Region.new("chr", 10, 15)
      reader.coverage(region).values.should eq([0, 0, 0, 0, 0])
      options = Mopdepth::Options.new(mode: Mopdepth::Mode::Fragment)
      reader.coverage(region, options: options).values.should eq([1, 1, 1, 1, 1])
      results = [] of Mopdepth::Coverage
      reader.each_coverage([region, Mopdepth::Region.new("chr", 23, 25)], options: options) { |result| results << result }
      results.map(&.values).should eq([[1, 1, 1, 1, 1], [0, 0]])
    end
  end
end
