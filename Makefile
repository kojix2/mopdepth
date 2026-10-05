SHELL := /bin/sh

.DEFAULT_GOAL := release

CRYSTAL ?= crystal
SHARDS ?= shards
CARGO ?= cargo
GIT ?= git

BIN_DIR ?= bin
BIN ?= $(BIN_DIR)/mopdepth
SOURCE := src/depth.cr
CRYSTAL_CACHE_DIR ?= $(CURDIR)/.build/crystal-cache

D4_VERSION ?= v0.3.11
D4_REPOSITORY ?= https://github.com/38/d4-format.git
D4_BUILD_ROOT ?= $(CURDIR)/.build/d4
D4_SOURCE_DIR := $(D4_BUILD_ROOT)/d4-format-$(D4_VERSION)
D4_CARGO_TARGET_DIR := $(D4_BUILD_ROOT)/cargo-target
D4_STATIC_DIR := $(D4_BUILD_ROOT)/lib
D4_STATIC_ARCHIVE := $(D4_STATIC_DIR)/libd4binding.a
D4_PATCH := $(CURDIR)/ci/d4-static.patch
D4_PATCH_STAMP := $(D4_SOURCE_DIR)/.mopdepth-static-patch
DEPS_STAMP := $(CURDIR)/.build/shards-release.stamp
DEV_DEPS_STAMP := $(CURDIR)/.build/shards-development.stamp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
MACOSX_DEPLOYMENT_TARGET ?= 11.0
D4_CARGO_ENV := MACOSX_DEPLOYMENT_TARGET=$(MACOSX_DEPLOYMENT_TARGET)
D4_SYSTEM_LIBS := -framework Security -framework CoreFoundation -framework SystemConfiguration -lresolv
else
D4_CARGO_ENV :=
D4_SYSTEM_LIBS := -ldl -lpthread -lm
endif

.PHONY: all release debug d4 deps dev-deps test test-d4 lint format format-check clean clean-d4

all: release

# Performance is a core requirement, so the default build is optimized.
release: deps
	@mkdir -p $(BIN_DIR) $(CRYSTAL_CACHE_DIR)
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) build $(SOURCE) --release -o $(BIN)

debug: deps
	@mkdir -p $(BIN_DIR) $(CRYSTAL_CACHE_DIR)
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) build $(SOURCE) -o $(BIN)

# Builds an optimized D4-enabled binary. libd4binding is copied into an
# archive-only directory so @[Link("d4binding")] cannot select a dylib/so.
d4: deps $(D4_STATIC_ARCHIVE)
	@mkdir -p $(BIN_DIR) $(CRYSTAL_CACHE_DIR)
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) build $(SOURCE) --release -Dd4 -o $(BIN) \
		--link-flags="-L$(D4_STATIC_DIR) $(D4_SYSTEM_LIBS)"

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

test-d4: d4 dev-deps
	CRYSTAL_CACHE_DIR=$(CRYSTAL_CACHE_DIR) $(CRYSTAL) spec spec/d4_output_spec.cr -Dd4 \
		--link-flags="-L$(D4_STATIC_DIR) $(D4_SYSTEM_LIBS)"

lint: dev-deps
	bin/ameba src spec

format:
	$(CRYSTAL) tool format src spec

format-check:
	$(CRYSTAL) tool format --check src spec

$(D4_SOURCE_DIR)/.git:
	@mkdir -p $(D4_BUILD_ROOT)
	$(GIT) clone --depth 1 --branch $(D4_VERSION) $(D4_REPOSITORY) $(D4_SOURCE_DIR)

$(D4_PATCH_STAMP): $(D4_SOURCE_DIR)/.git $(D4_PATCH)
	@cd $(D4_SOURCE_DIR) && \
		if $(GIT) apply --unidiff-zero --reverse --check $(D4_PATCH) >/dev/null 2>&1; then \
			:; \
		else \
			$(GIT) apply --unidiff-zero --check $(D4_PATCH) && \
				$(GIT) apply --unidiff-zero $(D4_PATCH); \
		fi
	@touch $@

$(D4_STATIC_ARCHIVE): $(D4_PATCH_STAMP) Makefile
	$(D4_CARGO_ENV) CARGO_TARGET_DIR=$(D4_CARGO_TARGET_DIR) $(CARGO) build \
		--manifest-path $(D4_SOURCE_DIR)/Cargo.toml --release -p d4binding
	@mkdir -p $(D4_STATIC_DIR)
	cp $(D4_CARGO_TARGET_DIR)/release/libd4binding.a $(D4_STATIC_ARCHIVE)

clean-d4:
	rm -rf $(D4_BUILD_ROOT)

clean:
	rm -f $(BIN)
	rm -rf $(CURDIR)/.build
