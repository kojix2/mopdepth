module Depth::Core
  record Options,
    mapq : Int32 = 0,
    min_frag_len : Int32 = -1,
    max_frag_len : Int32 = Int32::MAX,
    exclude_flag : UInt16 = 1796_u16,
    include_flag : UInt16 = 0_u16,
    fast_mode : Bool = false,
    fragment_mode : Bool = false,
    read_groups : Array(String) = [] of String
end
