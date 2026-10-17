require "hts"

# Small indexed input shared by CLI smoke tests. Keeping it generated makes
# the suite runnable without checking out mosdepth's large fixture submodule.
module MopdepthSpecFixture
  def self.create_bam(path : String) : Nil
    header = HTS::Bam::Header.parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n")
    HTS::Bam.open(path, "wb") do |bam|
      bam.write_header(header)
      bam.write(HTS::Bam::Record.new(
        header, "pair", 99_u16, "chr1", 2_i64, 60_u8, "6M", "AAAAAA",
        [30_u8] * 6, 0, 4_i64, 8_i64
      ))
      bam.write(HTS::Bam::Record.new(
        header, "pair", 147_u16, "chr1", 4_i64, 60_u8, "6M", "AAAAAA",
        [30_u8] * 6, 0, 2_i64, -8_i64
      ))
    end
    HTS::Bam.build_index(path, verbose: false)
  end
end
