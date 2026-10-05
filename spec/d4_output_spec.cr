require "./spec_helper"
require "../src/depth/config"
require "../src/depth/runner"

{% if flag?(:d4) %}
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
        file.chromosomes.should eq({"chr1" => 10_u32, "empty" => 5_u32})
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
{% else %}
  describe "D4 output" do
    it "reports that the default binary does not include D4" do
      config = Depth::Config.new
      config.prefix = "out"
      config.path = "input.bam"
      config.d4 = true

      expect_raises(ArgumentError, /rebuild with `make d4`/) do
        config.validate!
      end
    end
  end
{% end %}
