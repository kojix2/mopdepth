require "../core/target"
require "../errors"

{% if flag?(:d4) %}
  require "d4"
{% end %}

module Depth::FileIO
  {% if flag?(:d4) %}
    class D4Output
      DEFAULT_BUFFER_SIZE = 65_536

      getter path : String

      @writer : D4::Writer
      @intervals : Array(D4::Interval)
      @chromosome : String?
      @closed = false

      def initialize(@path : String, targets : Array(Core::Target), buffer_size : Int32 = DEFAULT_BUFFER_SIZE)
        raise Depth::OutputError.new("D4 interval buffer size must be positive") if buffer_size <= 0

        @writer = D4::Writer.new(@path)
        @intervals = Array(D4::Interval).new(buffer_size)
        @chromosome = nil
        @buffer_size = buffer_size
        begin
          chromosomes = targets.map { |target| {target.name, target.length.to_u32} }
          @writer.set_chromosomes(chromosomes)
        rescue ex
          @writer.close rescue nil
          raise ex
        end
      end

      def write_interval(chromosome : String, start : Int32, stop : Int32, value : Int32)
        raise Depth::OutputError.new("D4 output is closed") if @closed
        raise Depth::OutputError.new("Invalid D4 interval #{chromosome}:#{start}-#{stop}") if start < 0 || stop <= start

        if current = @chromosome
          flush if current != chromosome
        end
        @chromosome = chromosome
        @intervals << D4::Interval.new(start.to_u32, stop.to_u32, value)
        flush if @intervals.size >= @buffer_size
      end

      def close
        return if @closed

        failures = [] of String
        begin
          flush
        rescue ex
          failures << "flush failed: #{ex.message}"
        end
        begin
          @writer.close
        rescue ex
          failures << "close failed: #{ex.message}"
        end
        @closed = true
        raise Depth::OutputError.new(failures.join("; ")) unless failures.empty?
      end

      def build_index
        raise Depth::OutputError.new("D4 output must be closed before indexing") unless @closed

        D4.build_sfi_index(@path)
      end

      private def flush
        return if @intervals.empty?
        chromosome = @chromosome || raise Depth::OutputError.new("D4 interval chromosome is missing")

        @writer.write_intervals(chromosome, @intervals)
        @intervals.clear
      end
    end
  {% else %}
    class D4Output
      getter path : String

      def initialize(@path : String, targets : Array(Core::Target), buffer_size : Int32 = 0)
        raise Depth::OutputError.new("D4 output is not available in this build; rebuild with `make d4`")
      end

      def write_interval(chromosome : String, start : Int32, stop : Int32, value : Int32)
        raise Depth::OutputError.new("D4 output is not available in this build")
      end

      def close
      end

      def build_index
        raise Depth::OutputError.new("D4 output is not available in this build")
      end
    end
  {% end %}
end
