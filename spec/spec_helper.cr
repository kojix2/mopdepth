require "spec"
require "process"
require "file_utils"
require "compress/gzip"

# Helper to use the compiled binary instead of `crystal run` for speed
module TestBin
  @@built = false

  def self.binary : String
    ENV["MOPDEPTH_BIN"]? || File.expand_path("./bin/mopdepth")
  end

  def self.ensure_built!
    return if @@built && File.exists?(binary)
    unless File.exists?(binary)
      status = Process.run("shards", ["build", "--release"],
        chdir: Dir.current,
        output: Process::Redirect::Close,
        error: Process::Redirect::Close)
      raise "Failed to build with 'shards build'" unless status.success?
    end
    @@built = true
  end

  def self.run(args : Array(String), prefix : String, bam : String, temp_dir : String) : Process::Status
    ensure_built!
    Process.run(binary, args + [prefix, bam],
      chdir: Dir.current,
      env: {"PWD" => temp_dir},
      output: Process::Redirect::Close,
      # error: Process::Redirect::Close)
      error: Process::Redirect::Inherit)
  end
end

# IO helpers for reading text from plain or gzip files
module TestIO
  def self.read_text(path : String) : String
    if path.ends_with?(".gz")
      File.open(path) do |file|
        Compress::Gzip::Reader.open(file) do |gzip|
          gzip.gets_to_end
        end
      end
    else
      File.read(path)
    end
  end
end
