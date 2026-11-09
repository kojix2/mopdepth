module Mopdepth
  # Excludes a record only when a string-valued AUX tag matches exactly.
  record TagFilter, tag : String, value : String do
    def initialize(@tag : String, @value : String)
      unless tag.matches?(/\A[A-Za-z][A-Za-z0-9]\z/)
        raise ArgumentError.new("--exclude-tag expects TAG:VALUE with a SAM tag matching [A-Za-z][A-Za-z0-9]: #{tag.inspect}")
      end
      unless value.each_byte.all? { |byte| 32 <= byte <= 126 }
        raise ArgumentError.new("--exclude-tag values must contain printable ASCII characters")
      end
    end

    def self.parse_list(specification : String) : Array(self)
      specification.split(',').map do |item|
        tag, separator, value = item.partition(':')
        if separator.empty?
          raise ArgumentError.new("--exclude-tag expects TAG:VALUE: #{item.inspect}")
        end
        new(tag, value)
      end
    end
  end
end
