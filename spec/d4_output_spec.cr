require "./spec_helper"
require "../src/depth/config"
require "../src/depth/runner"
require "d4"

private def create_d4_test_bam(path : String)
  header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:empty\tLN:5\n")
  HTS::Bam.open(path, "wb") do |bam|
    bam.write_header(header)
    record = HTS::Bam::Record.new(
      header, "read", 0_u16, "chr1", 2_i64, 60_u8, "3M", "AAA",
      [30_u8, 30_u8, 30_u8], -1, 0_i64, 0_i64
    )
    bam.write(record)
  end
  HTS::Bam.build_index(path, verbose: false)
end

describe "D4 output" do
  temp_dir = "/tmp/mopdepth_d4_output"

  before_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "preserves depths at and above the dictionary boundary" do
    path = "#{temp_dir}/high-depth.d4"
    targets = [Depth::Core::Target.new("chr1", 5_i32, 0_i32)]
    output = Depth::FileIO::D4Output.new(path, targets)
    output.write_interval("chr1", 0_i32, 1_i32, 255_i32)
    output.write_interval("chr1", 1_i32, 2_i32, 256_i32)
    output.write_interval("chr1", 2_i32, 3_i32, 10_000_i32)
    output.close
    output.build_index

    D4.open(path) do |file|
      file.values("chr1").should eq([255, 256, 10_000, 0, 0])
    end
  end

  it "writes indexed per-base coverage without a BGZF per-base file" do
    bam_path = "#{temp_dir}/input.bam"
    create_d4_test_bam(bam_path)
    prefix = "#{temp_dir}/result"
    config = Depth::Config.new
    config.prefix = prefix
    config.path = bam_path
    config.d4 = true

    Depth::Runner.new(config).run

    File.exists?("#{prefix}.per-base.d4").should be_true
    File.exists?("#{prefix}.per-base.bed.gz").should be_false
    D4.open("#{prefix}.per-base.d4") do |file|
      file.chromosomes.map { |chromosome| {chromosome.name, chromosome.size} }.should eq([{"chr1", 10_i64}, {"empty", 5_i64}])
      file.default_track.has_index?(D4::IndexKind::SecondaryFrames).should be_true
      file.values("chr1").should eq([0, 0, 1, 1, 1, 0, 0, 0, 0, 0])
      file.values("empty").should eq([0, 0, 0, 0, 0])
    end
  end

  it "keeps positions outside a selected range at zero" do
    bam_path = "#{temp_dir}/partial.bam"
    create_d4_test_bam(bam_path)
    prefix = "#{temp_dir}/partial"
    config = Depth::Config.new
    config.prefix = prefix
    config.path = bam_path
    config.chrom = "chr1:3-4"
    config.d4 = true

    Depth::Runner.new(config).run

    D4.open("#{prefix}.per-base.d4") do |file|
      file.values("chr1").should eq([0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
      file.values("empty").should eq([0, 0, 0, 0, 0])
    end
  end
end
