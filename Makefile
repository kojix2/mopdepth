SHELL := /bin/sh

.DEFAULT_GOAL := release

CRYSTAL ?= crystal
SHARDS ?= shards
BIN_DIR ?= bin
BIN ?= $(BIN_DIR)/mopdepth
SOURCE := src/depth.cr
CRYSTAL_CACHE_DIR ?= $(CURDIR)/.build/crystal-cache

DEPS_STAMP := $(CURDIR)/.build/shards-release.stamp
DEV_DEPS_STAMP := $(CURDIR)/.build/shards-development.stamp

.PHONY: all release debug deps dev-deps test lint format format-check clean

all: release

# Performance is a core requirement, so the default build is optimized.
release: deps
	@mkdir -p $(BIN_DIR) $(CRYSTAL_CACHE_DIR)
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) build $(SOURCE) --release -o $(BIN)

debug: deps
	@mkdir -p $(BIN_DIR) $(CRYSTAL_CACHE_DIR)
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) build $(SOURCE) -o $(BIN)

deps: $(DEPS_STAMP)

dev-deps: $(DEV_DEPS_STAMP)

$(DEPS_STAMP): shard.yml shard.lock
	$(SHARDS) install --production
	@mkdir -p $(dir $@)
	@touch $@

$(DEV_DEPS_STAMP): shard.yml shard.lock
	$(SHARDS) install --frozen
	@mkdir -p $(dir $@)
	@touch $@

test: dev-deps
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) spec

lint: dev-deps
	bin/ameba src spec

format:
	$(CRYSTAL) tool format src spec

format-check:
	$(CRYSTAL) tool format --check src spec

clean:
	rm -f $(BIN)
	rm -rf $(CURDIR)/.build
