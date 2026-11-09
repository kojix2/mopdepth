require "./spec_helper"
require "../src/mopdepth"

describe "alignment read failures" do
  it "propagates a corrupt BAM data block and leaves outputs unindexed" do
    directory = File.tempname("mopdepth-corrupt-input")
    FileUtils.mkdir_p(directory)
    begin
      path = "#{directory}/input.bam"
      header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n")
      HTS::Bam.open(path, "wb") do |bam|
        bam.write_header(header)
        2_000.times do |index|
          record = HTS::Bam::Record.new(header, "read#{index}", 0_u16, "chr1", 2_i64,
            60_u8, "50M", "A" * 50, [30_u8] * 50, -1, 0_i64, 0_i64)
          bam.write(record)
        end
      end
      HTS::Bam.build_index(path, verbose: false)
      data = File.read(path).to_slice.dup
      # Leave the header and index intact, but invalidate a later BGZF CRC.
      second = data[16].to_i32 + (data[17].to_i32 << 8) + 1
      size = data[second + 16].to_i32 + (data[second + 17].to_i32 << 8) + 1
      data[second + size - 8] ^= 1_u8
      File.write(path, data)

      Mopdepth::Reader.open(path) do |reader|
        expect_raises(HTS::Bam::ReadError) { reader.coverage("chr1") }
      end
      TestBin.ensure_built!
      prefix = "#{directory}/result"
      error = IO::Memory.new
      status = Process.run(TestBin.binary, [prefix, path], output: Process::Redirect::Close, error: error)
      status.success?.should be_false
      error.to_s.should contain("Failed to read SAM/BAM")
      File.exists?("#{prefix}.per-base.bed.gz.csi").should be_false
      File.read("#{prefix}.mopdepth.summary.txt").should_not contain("total\t")
    ensure
      FileUtils.rm_rf(directory)
    end
  end
end
