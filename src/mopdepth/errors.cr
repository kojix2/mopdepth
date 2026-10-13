module Mopdepth
  class Error < Exception; end

  class InputError < Error; end

  class IndexError < Error; end

  class UnknownReferenceError < Error; end

  class ClosedError < Error; end

  class CancelledError < Error; end

  class CoverageOverflowError < Error; end
end
