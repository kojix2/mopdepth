module Depth::Stats
  # Parse threshold arguments string like "1,2,3" into array of integers
  def self.threshold_args(ts : String) : Array(Int32)
    return [] of Int32 if ts.empty? || ts == "nil"

    begin
      result = ts.split(',').map(&.to_i)
      result.sort!
      result
    rescue ArgumentError
      raise ArgumentError.new("Invalid threshold string: '#{ts}'")
    end
  end
end
