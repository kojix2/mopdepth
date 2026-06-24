module Depth::Core
  # coverage array diff representation and cumulative coverage alias
  alias Coverage = Array(Int32)

  struct DepthSegment
    getter start : Int32
    getter stop : Int32
    getter depth : Int32

    def initialize(@start : Int32, @stop : Int32, @depth : Int32)
    end

    def length : Int32
      @stop - @start
    end

    def empty? : Bool
      @stop <= @start
    end
  end
end
