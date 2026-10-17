require "./spec_helper"
require "../src/depth/config"
require "../src/depth/io/output_manager"
require "../src/depth/stats/distribution"

private def with_output_environment(mop_precision : String? = nil, mos_precision : String? = nil,
                                    bgzf_buffer : String? = nil, &)
  old_mop_precision = ENV["MOPDEPTH_PRECISION"]?
  old_mos_precision = ENV["MOSDEPTH_PRECISION"]?
  old_bgzf_buffer = ENV["MOPDEPTH_BGZF_BUFFER"]?

  set_output_environment("MOPDEPTH_PRECISION", mop_precision)
  set_output_environment("MOSDEPTH_PRECISION", mos_precision)
  set_output_environment("MOPDEPTH_BGZF_BUFFER", bgzf_buffer)
  yield
ensure
  set_output_environment("MOPDEPTH_PRECISION", old_mop_precision)
  set_output_environment("MOSDEPTH_PRECISION", old_mos_precision)
  set_output_environment("MOPDEPTH_BGZF_BUFFER", old_bgzf_buffer)
end

private def set_output_environment(name : String, value : String?)
  if value
    ENV[name] = value
  else
    ENV.delete(name)
  end
end

class OutputSettingsDistributionHost
  extend Depth::Stats::Distribution
end

describe "output environment settings" do
  temp_dir = File.tempname("mopdepth_output_settings")

  before_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  it "uses the same bounded precision for distribution output" do
    with_output_environment(mop_precision: "0") do
      io = IO::Memory.new
      OutputSettingsDistributionHost.write_distribution(io, "chr1", [0_i64, 1_i64])
      io.to_s.should eq("chr1\t1\t1\nchr1\t0\t1\n")
      Depth::OutputSettings.precision.should eq(0)
    end

    with_output_environment(mop_precision: Depth::OutputSettings::MAX_PRECISION.to_s) do
      Depth::OutputSettings.precision.should eq(Depth::OutputSettings::MAX_PRECISION)
    end

    ["-1", "19", "not-an-integer"].each do |value|
      with_output_environment(mop_precision: value) do
        expect_raises(Depth::ConfigError, /MOPDEPTH_PRECISION must be an integer between 0 and 18/) do
          OutputSettingsDistributionHost.write_distribution(IO::Memory.new, "chr1", [1_i64])
        end
      end
    end

    with_output_environment(mos_precision: "-1") do
      expect_raises(Depth::ConfigError, /MOSDEPTH_PRECISION must be an integer between 0 and 18/) do
        OutputSettingsDistributionHost.write_distribution(IO::Memory.new, "chr1", [1_i64])
      end
    end
  end

  it "writes a sparse distribution in descending order without changing it" do
    dist = Array(Int64).new(303, 0_i64)
    dist[302] = 1_i64
    original = dist.dup
    io = IO::Memory.new

    with_output_environment(mop_precision: "3") do
      OutputSettingsDistributionHost.write_distribution(io, "chr1", dist)
    end

    lines = io.to_s.lines
    lines.size.should eq(302)
    lines.first.should eq("chr1\t302\t1.000")
    lines[1].should eq("chr1\t300\t1.000")
    lines.last.should eq("chr1\t0\t1.000")
    dist.should eq(original)
  end

  it "rejects invalid precision before OutputManager creates output" do
    prefix = "#{temp_dir}/invalid-precision"
    config = Depth::Config.new
    config.prefix = prefix
    config.path = "unused.bam"

    with_output_environment(mop_precision: "-1") do
      expect_raises(Depth::ConfigError, /MOPDEPTH_PRECISION must be an integer between 0 and 18/) do
        config.validate!
      end
      expect_raises(Depth::ConfigError, /MOPDEPTH_PRECISION must be an integer between 0 and 18/) do
        Depth::FileIO::OutputManager.new(config)
      end
    end

    File.exists?("#{prefix}.mopdepth.summary.txt").should be_false
  end

  it "validates BGZF buffer limits and flushes with the minimum buffer" do
    with_output_environment(bgzf_buffer: "1") do
      Depth::OutputSettings.bgzf_buffer_size.should eq(1)
    end

    too_large = (Depth::OutputSettings::MAX_BGZF_BUFFER_SIZE + 1).to_s
    ["0", "-1", "not-an-integer", too_large].each do |value|
      with_output_environment(bgzf_buffer: value) do
        expect_raises(Depth::ConfigError, /MOPDEPTH_BGZF_BUFFER must be an integer between 1 and 16777216/) do
          Depth::OutputSettings.bgzf_buffer_size
        end
      end
    end

    prefix = "#{temp_dir}/tiny-buffer"
    config = Depth::Config.new
    config.prefix = prefix
    config.path = "unused.bam"
    config.by = "1"

    with_output_environment(bgzf_buffer: "1") do
      output = Depth::FileIO::OutputManager.new(config)
      output.write_per_base_interval("chr1", 0, 1, 7)
      output.close_all
    end

    TestIO.read_text("#{prefix}.per-base.bed.gz").should eq("chr1\t0\t1\t7\n")
  end
end
