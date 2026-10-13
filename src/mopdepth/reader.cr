require "hts"
require "../depth/core/options"
require "../depth/core/coverage"
require "../depth/core/region"
require "../depth/core/coverage_result"
require "../depth/core/coverage_calculator"

module Mopdepth
  # An indexed BAM/CRAM coverage reader.
  #
  # A Reader owns its HTS handle and is intended for serial use. Create one
  # Reader per concurrent worker. Call `close`, or use the block form of
  # `Reader.open`, to release the native handle deterministically.
  class Reader
    CHECK_INTERVAL = 65_536

    getter path : String
    getter targets : Array(Target)
    getter file_format : String
    getter records_examined = 0_i64

    @closed = false

    private def initialize(@path : String, @bam : HTS::Bam, @targets : Array(Target), @file_format : String)
      @targets_by_name = Hash(String, Target).new
      @targets.each { |target| @targets_by_name[target.name] = target }
    end

    def self.open(path : String | Path, *, reference : String | Path = "", threads : Int32 = 0) : self
      raise ArgumentError.new("threads must not be negative") if threads < 0
      input = path.to_s
      bam = nil.as(HTS::Bam?)
      begin
        bam = HTS::Bam.open(input, fai: reference.to_s, threads: threads)
        format = bam.file_format
        unless format == "bam" || format == "cram"
          raise InputError.new("indexed BAM or CRAM input is required; got #{format}")
        end
        begin
          bam.load_index
        rescue ex
          raise IndexError.new("failed to load the BAM/CRAM index for #{input}: #{ex.message}")
        end
        names = bam.header.target_names
        lengths = bam.header.target_len
        targets = names.map_with_index do |name, index|
          length = lengths[index].to_i64
          unless 0 < length < Int32::MAX
            raise InputError.new("reference length is outside the supported range for #{name}: #{length}")
          end
          Target.new(name, length.to_i32, index)
        end
        raise InputError.new("alignment header has no reference sequences") if targets.empty?
        new(input, bam, targets, format)
      rescue ex : Error
        bam.try(&.close)
        raise ex
      rescue ex
        bam.try(&.close)
        raise InputError.new("failed to open alignment #{input}: #{ex.message}")
      end
    end

    def self.open(path : String | Path, *, reference : String | Path = "", threads : Int32 = 0,
                  & : self -> T) : T forall T
      reader = open(path, reference: reference, threads: threads)
      begin
        yield reader
      ensure
        reader.close
      end
    end

    def target(name : String) : Target
      ensure_open
      @targets_by_name[name]? || raise UnknownReferenceError.new("reference not found: #{name}")
    end

    def coverage(chromosome : String, *, options : Options = Options.new,
                 control : Control = Control.new) : Coverage
      coverage(target(chromosome).region, options: options, control: control)
    end

    def coverage(region : Region, *, options : Options = Options.new,
                 control : Control = Control.new) : Coverage
      coverage_into(region, [] of Int32, options: options, control: control)
    end

    # Calculates dense coverage into a reusable caller-owned array. On return,
    # the array has exactly `region.length` values and contains no diff sentinel.
    def coverage_into(region : Region, values : Array(Int32), *, options : Options = Options.new,
                      control : Control = Control.new) : Coverage
      calculate(region, values, options, control, 1, 1)
    end

    def each_coverage(*, options : Options = Options.new, control : Control = Control.new,
                      & : Coverage ->) : Nil
      @targets.each_with_index do |item, index|
        yield calculate(item.region, [] of Int32, options, control, index + 1, @targets.size)
      end
    end

    def each_coverage(regions : Enumerable(Region), *, options : Options = Options.new,
                      control : Control = Control.new, & : Coverage ->) : Nil
      selected = regions.to_a
      selected.each_with_index do |region, index|
        yield calculate(region, [] of Int32, options, control, index + 1, selected.size)
      end
    end

    def close : Nil
      return if @closed
      @bam.close
      @closed = true
    end

    def closed? : Bool
      @closed
    end

    private def calculate(region : Region, values : Array(Int32), options : Options, control : Control,
                          target_index : Int32, target_count : Int32) : Coverage
      ensure_open
      options.validate!
      reference = target(region.chromosome)
      if region.stop > reference.length
        raise ArgumentError.new("region exceeds #{region.chromosome} length #{reference.length}: #{region.start}-#{region.stop}")
      end
      progress = Progress.new(Stage::Preparing, region, target_index, target_count, records_examined)
      control.update(progress)
      prepare_values(values, region.length + 1, control)
      return Coverage.new(region, values.tap(&.pop), CoverageStatus::NoData, 0_i64) if region.empty?

      core_options = Depth::Core::Options.new(
        mapq: options.minimum_mapping_quality,
        min_frag_len: options.minimum_fragment_length || -1,
        max_frag_len: options.maximum_fragment_length || Int32::MAX,
        exclude_flag: options.exclude_flags,
        include_flag: options.include_flags,
        fast_mode: options.mode.fast?,
        fragment_mode: options.mode.fragment?,
        read_groups: options.read_groups,
      )
      query = if options.mode.fragment?
                Depth::Core::Region.new(region.chromosome, 0, reference.length)
              else
                Depth::Core::Region.new(region.chromosome, region.start, region.stop)
              end
      before = records_examined
      progress = progress.copy_with(stage: Stage::Reading)
      control.update(progress)
      on_record = -> {
        @records_examined += 1
        control.update(progress.copy_with(records_examined: records_examined)) if records_examined % 8192 == 0
      }
      result = Depth::Core::CoverageCalculator.new(@bam, core_options, on_record).calculate(values, query, region.start)
      control.update(progress.copy_with(stage: Stage::Finalizing, records_examined: records_examined))
      status = result == Depth::Core::CoverageResult::NoData.value ? CoverageStatus::NoData : CoverageStatus::Observed
      finalize_values(values, region.length, status, control)
      Coverage.new(region, values, status, records_examined - before)
    end

    private def prepare_values(values : Array(Int32), size : Int32, control : Control) : Nil
      values.clear
      while values.size < size
        control.check!
        values.concat(Array(Int32).new(Math.min(CHECK_INTERVAL, size - values.size), 0_i32))
      end
    end

    private def finalize_values(values : Array(Int32), length : Int32, status : CoverageStatus,
                                control : Control) : Nil
      if status.observed?
        sum = 0_i64
        length.times do |index|
          control.check! if index % CHECK_INTERVAL == 0
          sum += values[index]
          unless 0 <= sum <= Int32::MAX
            raise CoverageOverflowError.new("coverage exceeds Int32 at offset #{index}: #{sum}")
          end
          values[index] = sum.to_i32
        end
      end
      values.pop
    end

    private def ensure_open : Nil
      raise ClosedError.new("mopdepth reader is closed") if closed?
    end
  end
end
