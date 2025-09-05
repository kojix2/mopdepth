module Depth::Core
  struct Region
    getter chrom : String
    getter start : Int32
    getter stop : Int32
    getter name : String?

    def initialize(@chrom : String, @start : Int32, @stop : Int32, @name : String? = nil)
      raise ArgumentError.new("start > stop for #{chrom}:#{start}-#{stop}") if @start > @stop
    end

    def to_s(io)
      if name
        io << chrom << ":" << start + 1 << "-" << stop << "\t" << name
      else
        io << chrom << ":" << start + 1 << "-" << stop
      end
    end
  end
end
