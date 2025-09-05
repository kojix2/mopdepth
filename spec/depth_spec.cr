require "./spec_helper"
require "../src/depth/version"

describe Depth do
  it "has a version number" do
    Depth::VERSION.should be_a(String)
  end
end
