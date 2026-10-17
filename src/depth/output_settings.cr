require "./errors"

module Depth
  # Environment-controlled output settings shared by every output path.
  module OutputSettings
    DEFAULT_PRECISION        =         2_i32
    MAX_PRECISION            =        18_i32
    DEFAULT_BGZF_BUFFER_SIZE = 2_097_152_i32
    MAX_BGZF_BUFFER_SIZE     = 16_i32 * 1024 * 1024

    def self.precision : Int32
      integer_env(
        "MOPDEPTH_PRECISION",
        fallback: "MOSDEPTH_PRECISION",
        default: DEFAULT_PRECISION,
        min: 0_i32,
        max: MAX_PRECISION,
      )
    end

    def self.bgzf_buffer_size : Int32
      integer_env(
        "MOPDEPTH_BGZF_BUFFER",
        default: DEFAULT_BGZF_BUFFER_SIZE,
        min: 1_i32,
        max: MAX_BGZF_BUFFER_SIZE,
      )
    end

    def self.validate! : Nil
      precision
      bgzf_buffer_size
    end

    private def self.integer_env(name : String, *, fallback : String? = nil,
                                 default : Int32, min : Int32, max : Int32) : Int32
      selected_name = name
      value = ENV[name]?
      if value.nil? && fallback
        selected_name = fallback
        value = ENV[fallback]?
      end
      return default unless value

      parsed = value.to_i32?
      unless parsed && min <= parsed <= max
        raise ConfigError.new("#{selected_name} must be an integer between #{min} and #{max}; got #{value.inspect}")
      end
      parsed
    end
  end
end
