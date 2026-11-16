MAGMA ?= magma
MAGMA_PATH := $(shell command -v "$(MAGMA)")

CHECK_FILES := $(patsubst tests/%,check/%,$(filter-out tests/_%,$(wildcard tests/*.m))) \
               $(patsubst examples/%,check/%,$(filter-out examples/_%,$(wildcard examples/*.m)))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

check/%.m: tests/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$< run_tests.m
check/%.m: examples/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$< run_tests.m

debug/%.m: tests/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$< debug:='' run_tests.m
debug/%.m: examples/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$< debug:='' run_tests.m

check: $(CHECK_FILES)

check-no-dups:
	@sh scripts/no_dup_test_basenames.sh

.PHONY: all check check-no-dups
