MAGMA ?= magma
MAGMA_PATH := $(shell command -v "$(MAGMA)")

CHECK_FILES := $(patsubst Tests/%,check/%,$(filter-out Tests/_%,$(wildcard Tests/*.m))) \
               $(patsubst examples/%,check/%,$(filter-out examples/_%,$(wildcard examples/*.m)))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

check/%.m: Tests/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -b filename:=$< exitsignal:='' run_tests.m
check/%.m: examples/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -b filename:=$< exitsignal:='' run_tests.m

debug/%.m: Tests/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -b filename:=$< debug:='' run_tests.m
debug/%.m: examples/%.m | check-no-dups
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -b filename:=$< debug:='' run_tests.m

check: $(CHECK_FILES)

check-no-dups:
	@sh scripts/no_dup_test_basenames.sh

.PHONY: all check check-no-dups

REFERENCE_PROCESSES ?= 8

REFERENCE_FILES := $(wildcard reference/checks/*.m reference/f3/checks/*.m)
HEAVY_FILES := $(wildcard reference/heavy/check_*.m)

export HEAVY_PREFIX_INDICES HEAVY_TARGET_LIMIT

check-heavy:
	@set -e; for check in $(HEAVY_FILES); do $(MAKE) reference-check/$$check; done

check-reference:
	@set -e; for check in $(REFERENCE_FILES); do $(MAKE) reference-check/$$check; done

reference-check/%.m:
	@mkdir -p artifacts/test-tmp
	@result=$$(mktemp "$(CURDIR)/artifacts/reference-result-XXXXXXXX"); \
	log=$$(mktemp "$(CURDIR)/artifacts/reference-log-XXXXXXXX"); \
	trap 'rm -f "$$result" "$$log"' EXIT HUP INT TERM; \
	check="$*"; case "$$check" in reference/*) ;; *) \
	  if test -f "reference/checks/$$check.m"; then check="reference/checks/$$check"; \
	  else check="reference/f3/checks/$$check"; fi ;; esac; \
	REFERENCE_PROCESSES="$(REFERENCE_PROCESSES)" MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" sh scripts/run_reference_check.sh "$$check.m" "$$result" "$$log" || { cat "$$log"; exit 1; }

check-slow:
	@set -e; for check in examples/slow/*.m; do \
	  mkdir -p artifacts/test-tmp; \
	  MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -b filename:=$$check exitsignal:='' run_tests.m; \
	done

.PHONY: check-reference check-slow check-heavy
