require "./spec_helper"
require "../src/depth/config"
require "../src/mopdepth"

describe Mopdepth::TagFilter do
  it "parses exact string matches, including empty values and colons" do
    filters = Mopdepth::TagFilter.parse_list("DT:SQ,X1:a:b,ZZ:")
    filters.map { |filter| {filter.tag, filter.value} }.should eq([
      {"DT", "SQ"}, {"X1", "a:b"}, {"ZZ", ""},
    ])
  end

  it "rejects malformed lists and invalid SAM tag names" do
    ["", "DT", "D:SQ", "DTT:SQ", "1T:SQ", "D_:SQ", "éT:SQ",
     "DT:SQ,", ",DT:SQ", "DT:SQ,,DT:LB"].each do |specification|
      expect_raises(ArgumentError, /--exclude-tag/) { Mopdepth::TagFilter.parse_list(specification) }
    end
    ["a\tb", "a\nb", "a\0b", "日本語"].each do |value|
      expect_raises(ArgumentError) { Mopdepth::TagFilter.new("DT", value) }
    end
  end

  it "passes the same validated filters to CLI and library options" do
    filters = Mopdepth::TagFilter.parse_list("DT:SQ")
    config = Depth::Config.new
    config.exclude_tags = filters
    config.to_options.exclude_tags.should eq(filters)
    Mopdepth::Options.new(exclude_tags: filters).exclude_tags.should eq(filters)
  end
end
