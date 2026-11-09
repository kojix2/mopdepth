require "./spec_helper"
require "../src/mopdepth"
require "d4"

private def create_tagged_alignment(path : String, mode : String = "wb", reference : String = "")
  header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:empty\tLN:5\n@RG\tID:keep\tSM:test\n@RG\tID:drop\tSM:test\n")
  HTS::Bam.open(path, mode, fai: reference) do |bam|
    bam.write_header(header)
    [nil, "SQ", "LB"].each_with_index do |tag, index|
      flag = index == 0 ? 0_u16 : 1024_u16
      record = HTS::Bam::Record.new(header, "read#{index}", flag, "chr1", 2_i64,
        index == 2 ? 10_u8 : 60_u8, "3M", "AAA", [30_u8] * 3, -1, 0_i64, 0_i64)
      record.aux.update_string("DT", tag) if tag
      record.aux.update_string("RG", index == 2 ? "drop" : "keep")
      bam.write(record)
    end
  end
  HTS::Bam.build_index(path, verbose: false)
end

describe "tag exclusion" do
  temp_dir = "/tmp/mopdepth_tag_exclusion"

  before_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "combines flag filters with absent, multiple, and repeated tag conditions" do
    bam_path = "#{temp_dir}/input.bam"
    create_tagged_alignment(bam_path)
    cases = [
      {[] of String, 1},
      {["-F", "772"], 3},
      {["-F", "772", "--exclude-tag", "DT:SQ"], 2},
      {["-F", "772", "--exclude-tag", "ZZ:nope"], 3},
      {["-F", "772", "--exclude-tag", "DT:SQ,DT:LB"], 1},
      {["-F", "772", "--exclude-tag", "DT:SQ", "--exclude-tag", "DT:LB"], 1},
      {["-F", "772", "--exclude-tag", "DT:SQ", "-Q", "20"], 1},
      {["-F", "772", "--exclude-tag", "DT:SQ", "-R", "keep"], 1},
      {["-F", "772", "--exclude-tag", "DT:SQ", "-x"], 2},
    ]
    cases.each_with_index do |(arguments, depth), index|
      prefix = "#{temp_dir}/case#{index}"
      TestBin.run(arguments, prefix, bam_path, temp_dir).success?.should be_true
      TestIO.read_text("#{prefix}.per-base.bed.gz").should eq(
        "chr1\t0\t2\t0\nchr1\t2\t5\t#{depth}\nchr1\t5\t20\t0\nempty\t0\t5\t0\n")
    end
  end

  it "filters library coverage in standard, fast, and fragment modes" do
    path = "#{temp_dir}/pairs.bam"
    header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n")
    HTS::Bam.open(path, "wb") do |bam|
      bam.write_header(header)
      [{99_u16, 2_i64, 6_i64, 7_i64}, {147_u16, 6_i64, 2_i64, -7_i64}].each do |flag, position, mate, length|
        ["SQ", "LB"].each do |tag|
          record = HTS::Bam::Record.new(header, tag, flag, "chr1", position, 60_u8,
            "3M", "AAA", [30_u8] * 3, 0, mate, length)
          record.aux.update_string("DT", tag)
          bam.write(record)
        end
      end
    end
    HTS::Bam.build_index(path, verbose: false)
    Mopdepth::Reader.open(path) do |reader|
      Mopdepth::Mode.each do |mode|
        options = Mopdepth::Options.new(mode: mode, exclude_tags: Mopdepth::TagFilter.parse_list("DT:SQ"))
        expected = Array(Int32).new(20, 0)
        (2...5).each { |index| expected[index] = 1 }
        (6...9).each { |index| expected[index] = 1 }
        expected[5] = 1 if mode.fragment?
        reader.coverage("chr1", options: options).values.should eq(expected)
      end
    end
    prefix = "#{temp_dir}/fragment"
    TestBin.run(["-a", "--exclude-tag", "DT:SQ"], prefix, path, temp_dir).success?.should be_true
    TestIO.read_text("#{prefix}.per-base.bed.gz").should contain("chr1\t2\t9\t1\n")
  end

  it "matches exact string types without coercing numeric tags" do
    path = "#{temp_dir}/types.bam"
    header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n")
    HTS::Bam.open(path, "wb") do |bam|
      bam.write_header(header)
      ["7", "", "a:b", "SQ"].each_with_index do |value, index|
        record = HTS::Bam::Record.new(header, "read#{index}", 0_u16, "chr1", 2_i64,
          60_u8, "3M", "AAA", [30_u8] * 3, -1, 0_i64, 0_i64)
        record.aux.update_string("DT", value)
        record.aux.update_int("NM", 7)
        record.aux.update_char("XC", '7')
        record.aux.update_array("XA", [7])
        record.aux.update_hex("HX", "ABCD")
        bam.write(record)
      end
    end
    HTS::Bam.build_index(path, verbose: false)
    Mopdepth::Reader.open(path) do |reader|
      {"NM:7" => 4, "XC:7" => 4, "XA:7" => 4, "HX:ABCD" => 0,
       "DT:" => 3, "DT:a:b" => 3, "DT:sq" => 4,
       "DT:S" => 4, "DT:SQx" => 4}.each do |condition, depth|
        options = Mopdepth::Options.new(exclude_tags: Mopdepth::TagFilter.parse_list(condition))
        reader.coverage("chr1", options: options).values[2, 3].should eq([depth] * 3)
      end
    end
  end

  it "decodes tag filters from CRAM in the CLI and library" do
    fasta = "#{temp_dir}/reference.fa"
    File.write(fasta, ">chr1\n#{"A" * 20}\n>empty\nAAAAA\n")
    HTS::Faidx.build_index(fasta)
    path = "#{temp_dir}/input.cram"
    create_tagged_alignment(path, "wc", fasta)
    options = Mopdepth::Options.new(exclude_flags: 772_u16, exclude_tags: Mopdepth::TagFilter.parse_list("DT:SQ"))
    Mopdepth::Reader.open(path, reference: fasta) do |reader|
      reader.coverage("chr1", options: options).values[2, 3].should eq([2, 2, 2])
    end
    prefix = "#{temp_dir}/cram"
    TestBin.run(["-f", fasta, "-F", "772", "--exclude-tag", "DT:SQ"], prefix, path, temp_dir).success?.should be_true
    TestIO.read_text("#{prefix}.per-base.bed.gz").should contain("chr1\t2\t5\t2\n")
  end

  it "applies exclusion to summaries, D4, quantization, and thresholds" do
    path = "#{temp_dir}/input.bam"
    create_tagged_alignment(path)
    prefix = "#{temp_dir}/outputs"
    TestBin.run(["-F", "772", "--exclude-tag", "DT:SQ", "--d4", "-q", "0:1:3:", "-b", "10", "-T", "0,1,2"],
      prefix, path, temp_dir).success?.should be_true
    File.read("#{prefix}.mopdepth.summary.txt").should contain("chr1\t20\t6\t0.30\t0\t2\n")
    TestIO.read_text("#{prefix}.quantized.bed.gz").should contain("chr1\t2\t5\t1:3\n")
    TestIO.read_text("#{prefix}.thresholds.bed.gz").should contain("chr1\t0\t10\tunknown\t10\t3\t3\n")
    D4.open("#{prefix}.per-base.d4") { |file| file.values("chr1")[2, 3].should eq([2, 2, 2]) }
  end

  it "rejects malformed conditions before opening the input or creating outputs" do
    TestBin.ensure_built!
    prefix = "#{temp_dir}/invalid"
    ["DT", "1T:SQ", "DT:SQ,"].each do |condition|
      error = IO::Memory.new
      status = Process.run(TestBin.binary, ["--exclude-tag", condition, prefix, "missing.bam"],
        output: Process::Redirect::Close, error: error)
      status.success?.should be_false
      error.to_s.should contain("--exclude-tag")
    end
    Dir.glob("#{prefix}*").should be_empty
  end
end
