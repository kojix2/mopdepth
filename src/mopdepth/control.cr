module Mopdepth
  # Cooperative cancellation and progress reporting for coverage operations.
  class Control
    @cancelled = Atomic(Bool).new(false)

    def initialize(@on_progress : Proc(Progress, Nil)? = nil)
    end

    def initialize(&@on_progress : Progress -> Nil)
    end

    def cancel : Nil
      @cancelled.set(true)
    end

    def cancelled? : Bool
      @cancelled.get
    end

    def check! : Nil
      raise CancelledError.new("coverage calculation was cancelled") if cancelled?
    end

    def update(progress : Progress) : Nil
      check!
      @on_progress.try(&.call(progress))
      check!
    end
  end
end
