require "./spec_helper"
require "../src/depth/config"
require "../src/depth/io/bed_reader"
require "../src/depth/io/output_manager"
require "../src/depth/stats/distribution"
require "hts"

private def create_audit_bam(path : String, header_text : String, & : HTS::Bam, HTS::Bam::Header ->)
  header = HTS::Bam::Header.parse(header_text)
  HTS::Bam.open(path, "wb") do |bam|
    bam.write_header(header)
    yield bam, header
  end
  HTS::Bam.build_index(path, verbose: false)
end

private def add_audit_record(bam : HTS::Bam, header : HTS::Bam::Header,
                             qname : String, chrom : String, pos : Int64, cigar : String, length : Int32,
                             mapq : UInt8 = 60_u8, flag : UInt16 = 0_u16,
                             mtid : Int32 = -1, mpos : Int64 = 0_i64, isize : Int64 = 0_i64)
  sequence = "A" * length
  record = HTS::Bam::Record.new(
    header, qname, flag, chrom, pos, mapq, cigar, sequence,
    Array(UInt8).new(length, 30_u8), mtid, mpos, isize
  )
  bam.write(record)
end

private def run_audit_mopdepth(args : Array(String), prefix : String, bam : String, temp_dir : String)
  status = TestBin.run(args, prefix, bam, temp_dir)
  status.success?.should be_true
end

class FailingIndexOutputManager < Depth::FileIO::OutputManager
  getter attempts = 0

  protected def build_csi_index(gz_path : String, csi : String) : Int32
    @attempts += 1
    -7
  end
end

class AuditDistributionHost
  extend Depth::Stats::Distribution
end

describe "audit regression coverage" do
  temp_dir = "/tmp/mopdepth_audit_regression"

  before_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "clears a reused coverage buffer across a short NoData target" do
    bam_path = "#{temp_dir}/reuse.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:A\tLN:100\n@SQ\tSN:B\tLN:10\n@SQ\tSN:C\tLN:100\n") do |bam, header|
      add_audit_record(bam, header, "a", "A", 20, "5M", 5)
      add_audit_record(bam, header, "c", "C", 0, "1M", 1)
    end

    prefix = "#{temp_dir}/reuse"
    run_audit_mopdepth([] of String, prefix, bam_path, temp_dir)
    output = TestIO.read_text("#{prefix}.per-base.bed.gz")
    output.should contain("C\t0\t1\t1\nC\t1\t100\t0\n")
    File.read("#{prefix}.mopdepth.summary.txt").should contain("C\t100\t1\t0.01\t0\t1")
  end

  it "applies template-length filtering in normal and fast modes" do
    bam_path = "#{temp_dir}/length.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n") do |bam, header|
      add_audit_record(bam, header, "pair", "chr1", 0, "10M", 10,
        flag: 99_u16, mtid: 0, mpos: 10, isize: 20)
      add_audit_record(bam, header, "pair", "chr1", 10, "10M", 10,
        flag: 147_u16, mtid: 0, mpos: 0, isize: -20)
    end

    [([] of String), ["-x"], ["-a"]].each_with_index do |mode, index|
      prefix = "#{temp_dir}/length#{index}"
      run_audit_mopdepth(mode + ["-l", "30"], prefix, bam_path, temp_dir)
      File.read("#{prefix}.mopdepth.summary.txt").should contain("chr1\t100\t0\t0.00\t0\t0")
    end
  end

  it "distinguishes an observed but fully filtered target from NoData" do
    bam_path = "#{temp_dir}/filtered.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:A\tLN:100\n@SQ\tSN:B\tLN:100\n") do |bam, header|
      add_audit_record(bam, header, "a", "A", 0, "10M", 10)
      add_audit_record(bam, header, "b", "B", 0, "10M", 10, mapq: 0_u8)
    end

    prefix = "#{temp_dir}/filtered"
    run_audit_mopdepth(["-Q", "20", "-b", "100", "-T", "0"], prefix, bam_path, temp_dir)
    File.read("#{prefix}.mopdepth.summary.txt").should contain("B\t100\t0\t0.00\t0\t0")
    TestIO.read_text("#{prefix}.thresholds.bed.gz").should contain("B\t0\t100\tunknown\t100")
  end

  it "does not skip BED-external targets when quantization needs them" do
    bam_path = "#{temp_dir}/skip.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:A\tLN:10\n@SQ\tSN:B\tLN:10\n") do |bam, header|
      add_audit_record(bam, header, "a", "A", 0, "1M", 1)
      add_audit_record(bam, header, "b", "B", 0, "1M", 1)
    end
    bed = "#{temp_dir}/only-a.bed"
    File.write(bed, "A\t0\t10\n")

    prefix = "#{temp_dir}/skip"
    run_audit_mopdepth(["-n", "-b", bed, "-q", "0:1:"], prefix, bam_path, temp_dir)
    TestIO.read_text("#{prefix}.quantized.bed.gz").should contain("B\t0\t1\t1:inf")
    File.read("#{prefix}.mopdepth.summary.txt").should contain("B\t10\t1")

    threshold_prefix = "#{temp_dir}/skip-threshold"
    run_audit_mopdepth(["-n", "-b", bed, "-T", "0"], threshold_prefix, bam_path, temp_dir)
    File.read("#{threshold_prefix}.mopdepth.summary.txt").should contain("B\t10\t1")
  end

  it "counts fragments whose reads lie outside a selected region" do
    bam_path = "#{temp_dir}/fragment.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n") do |bam, header|
      add_audit_record(bam, header, "pair", "chr1", 100, "100M", 100,
        flag: 99_u16, mtid: 0, mpos: 400, isize: 400)
      add_audit_record(bam, header, "pair", "chr1", 400, "100M", 100,
        flag: 147_u16, mtid: 0, mpos: 100, isize: -400)
    end

    prefix = "#{temp_dir}/fragment"
    run_audit_mopdepth(["-a", "-c", "chr1:251-350"], prefix, bam_path, temp_dir)
    TestIO.read_text("#{prefix}.per-base.bed.gz").should eq("chr1\t250\t350\t1\n")
  end

  it "keeps the union depth for a contained read pair" do
    bam_path = "#{temp_dir}/contained.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n") do |bam, header|
      add_audit_record(bam, header, "pair", "chr1", 100, "100M", 100,
        flag: 99_u16, mtid: 0, mpos: 120, isize: 100)
      add_audit_record(bam, header, "pair", "chr1", 120, "50M", 50,
        flag: 147_u16, mtid: 0, mpos: 100, isize: -100)
    end

    prefix = "#{temp_dir}/contained"
    run_audit_mopdepth([] of String, prefix, bam_path, temp_dir)
    File.read("#{prefix}.mopdepth.summary.txt").should contain("chr1\t300\t100")
    TestIO.read_text("#{prefix}.per-base.bed.gz").should contain("chr1\t100\t200\t1")
  end

  it "does not treat a single deletion CIGAR as covered sequence" do
    bam_path = "#{temp_dir}/deletion.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n") do |bam, header|
      add_audit_record(bam, header, "pair", "chr1", 100, "50D", 0,
        flag: 99_u16, mtid: 0, mpos: 120, isize: 70)
      add_audit_record(bam, header, "pair", "chr1", 120, "50M", 50,
        flag: 147_u16, mtid: 0, mpos: 100, isize: -70)
    end

    prefix = "#{temp_dir}/deletion"
    run_audit_mopdepth([] of String, prefix, bam_path, temp_dir)
    File.read("#{prefix}.mopdepth.summary.txt").should contain("chr1\t300\t50")
    TestIO.read_text("#{prefix}.per-base.bed.gz").should contain("chr1\t120\t170\t1")
  end

  it "handles equivalent split CIGARs and a mate filtered on only one side" do
    bam_path = "#{temp_dir}/mate-edges.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:equiv\tLN:300\n@SQ\tSN:filtered\tLN:100\n") do |bam, header|
      add_audit_record(bam, header, "equiv", "equiv", 100, "99M1=", 100,
        flag: 99_u16, mtid: 0, mpos: 120, isize: 100)
      add_audit_record(bam, header, "equiv", "equiv", 120, "50M", 50,
        flag: 147_u16, mtid: 0, mpos: 100, isize: -100)
      add_audit_record(bam, header, "filtered", "filtered", 10, "10M", 10,
        flag: 99_u16, mtid: 1, mpos: 15, isize: 15)
      add_audit_record(bam, header, "filtered", "filtered", 15, "10M", 10,
        mapq: 0_u8, flag: 147_u16, mtid: 1, mpos: 10, isize: -15)
    end

    prefix = "#{temp_dir}/mate-edges"
    run_audit_mopdepth(["-Q", "20"], prefix, bam_path, temp_dir)
    summary = File.read("#{prefix}.mopdepth.summary.txt")
    summary.should contain("equiv\t300\t100")
    summary.should contain("filtered\t100\t10")
  end

  it "rounds exact half-window means away from zero" do
    bam_path = "#{temp_dir}/rounding.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:half\tLN:2\n@SQ\tSN:twohalf\tLN:2\n") do |bam, header|
      add_audit_record(bam, header, "half", "half", 0, "1M", 1)
      2.times do |i|
        add_audit_record(bam, header, "base#{i}", "twohalf", 0, "2M", 2)
      end
      add_audit_record(bam, header, "extra", "twohalf", 0, "1M", 1)
    end

    prefix = "#{temp_dir}/rounding"
    run_audit_mopdepth(["-b", "2"], prefix, bam_path, temp_dir)
    dist = File.read("#{prefix}.mopdepth.region.dist.txt")
    dist.should contain("half\t1\t1.00")
    dist.should contain("twohalf\t3\t1.00")
  end

  it "always writes the summary header before empty totals" do
    bam_path = "#{temp_dir}/empty.bam"
    create_audit_bam(bam_path, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:empty\tLN:10\n") { |_bam, _header| }

    prefix = "#{temp_dir}/empty"
    run_audit_mopdepth([] of String, prefix, bam_path, temp_dir)
    File.read("#{prefix}.mopdepth.summary.txt").should eq(
      "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t0\t0\t0.00\t0\t0\n"
    )
  end

  it "reads gzip BED, ignores blank lines, and normalizes empty names" do
    bed_path = "#{temp_dir}/regions.bed.gz"
    File.open(bed_path, "w") do |file|
      Compress::Gzip::Writer.open(file) do |gzip|
        gzip << "chr1\t0\t2\t\n\nchr1\t4\t6\tname\n"
      end
    end
    targets = [Depth::Core::Target.new("chr1", 10, 0)]

    regions = Depth::FileIO.read_bed(bed_path, targets)["chr1"]
    regions.size.should eq(2)
    regions[0].name.should be_nil
    regions[1].name.should eq("name")
  end

  it "uses an explicit FASTA reference to decode CRAM" do
    fasta = "#{temp_dir}/reference.fa"
    File.write(fasta, ">chr1\n#{"A" * 100}\n")
    HTS::Faidx.build_index(fasta)

    cram_path = "#{temp_dir}/input.cram"
    header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n")
    HTS::Bam.open(cram_path, "wc", fai: fasta) do |cram|
      cram.write_header(header)
      add_audit_record(cram, header, "read", "chr1", 10, "5M", 5)
    end
    HTS::Bam.build_index(cram_path, verbose: false)

    prefix = "#{temp_dir}/cram"
    run_audit_mopdepth(["-f", fasta], prefix, cram_path, temp_dir)
    TestIO.read_text("#{prefix}.per-base.bed.gz").should contain("chr1\t10\t15\t1")

    env_prefix = "#{temp_dir}/cram-env"
    status = Process.run(TestBin.binary, [env_prefix, cram_path],
      env: {"REF_PATH" => fasta}, output: Process::Redirect::Close, error: Process::Redirect::Close)
    status.success?.should be_true
    TestIO.read_text("#{env_prefix}.per-base.bed.gz").should contain("chr1\t10\t15\t1")
  end

  it "rejects unsafe CLI arguments before opening input" do
    TestBin.ensure_built!
    [
      ["-b", "0", "prefix", "missing.bam"],
      ["-Q", "not-a-number", "prefix", "missing.bam"],
      ["prefix", "missing.bam", "extra"],
    ].each do |arguments|
      status = Process.run(TestBin.binary, arguments,
        output: Process::Redirect::Close, error: Process::Redirect::Close)
      status.success?.should be_false
    end
  end

  it "rejects invalid BED intervals and ignores unknown chromosomes" do
    targets = [Depth::Core::Target.new("chr1", 10, 0)]
    invalid = "#{temp_dir}/invalid.bed"
    File.write(invalid, "chr1\t0\t11\n")
    expect_raises(Depth::ConfigError, /invalid\.bed:1/) do
      Depth::FileIO.read_bed(invalid, targets)
    end

    unknown = "#{temp_dir}/unknown.bed"
    File.write(unknown, "missing\t0\t1\nmissing\t2\t3\nchr1\t0\t1\n")
    Depth::FileIO.read_bed(unknown, targets)["chr1"].size.should eq(1)
  end

  it "tries every CSI build and reports failures" do
    config = Depth::Config.new
    config.prefix = "#{temp_dir}/failure"
    config.path = "unused.bam"
    config.by = "10"
    output = FailingIndexOutputManager.new(config)

    error = expect_raises(Depth::OutputError) { output.close_all }
    output.attempts.should eq(2)
    error.message.to_s.should contain("ret=-7")
  end

  it "keeps the final real base in quantized output" do
    quants = Depth::Stats::Quantize.get_quantize_args("0:1:")
    coverage = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    segments = [] of Tuple(Int32, Int32, String)
    Depth::Stats::Quantize.gen_quantized(quants, coverage) { |segment| segments << segment }
    segments.should eq([
      {0, 9, "0:1"},
      {9, 10, "1:inf"},
    ])
  end

  it "treats an explicit Int32 maximum quantize boundary as finite" do
    quants = Depth::Stats::Quantize.get_quantize_args("0:#{Int32::MAX}:")
    quants.should eq([0_i64, Int32::MAX.to_i64, Int64::MAX])
    Depth::Stats::Quantize.make_lookup(quants).should eq([
      "0:#{Int32::MAX}",
      "#{Int32::MAX}:inf",
    ])
  end

  it "uses MOPDEPTH_PRECISION consistently for distributions" do
    old_mop = ENV["MOPDEPTH_PRECISION"]?
    old_mos = ENV["MOSDEPTH_PRECISION"]?
    ENV["MOPDEPTH_PRECISION"] = "3"
    ENV["MOSDEPTH_PRECISION"] = "1"
    io = IO::Memory.new
    AuditDistributionHost.write_distribution(io, "chr1", [0_i64, 1_i64])
    io.to_s.should eq("chr1\t1\t1.000\nchr1\t0\t1.000\n")
  ensure
    if old_mop
      ENV["MOPDEPTH_PRECISION"] = old_mop
    else
      ENV.delete("MOPDEPTH_PRECISION")
    end
    if old_mos
      ENV["MOSDEPTH_PRECISION"] = old_mos
    else
      ENV.delete("MOSDEPTH_PRECISION")
    end
  end

  it "falls back to the default quantize label when an environment label is empty" do
    old_label = ENV["MOSDEPTH_Q0"]?
    ENV["MOSDEPTH_Q0"] = ""
    Depth::Stats::Quantize.make_lookup([0_i64, 1_i64]).should eq(["0:1"])
  ensure
    if old_label
      ENV["MOSDEPTH_Q0"] = old_label
    else
      ENV.delete("MOSDEPTH_Q0")
    end
  end
end
