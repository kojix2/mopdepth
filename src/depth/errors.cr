module Depth
  class Error < Exception; end

  class FileNotFoundError < Error; end

  class InvalidRegionError < Error; end

  class CoverageCalculationError < Error; end

  class ConfigError < Error; end

  class BamIndexError < Error; end
end
