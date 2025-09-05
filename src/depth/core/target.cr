module Depth::Core
  struct Target
    getter name : String
    getter length : Int32
    getter tid : Int32

    def initialize(@name : String, @length : Int32, @tid : Int32)
    end
  end
end
