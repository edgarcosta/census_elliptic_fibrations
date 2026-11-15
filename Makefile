MAGMA ?= magma
MAGMA_PATH := $(shell command -v "$(MAGMA)")

CXXFLAGS ?= -O2 -Wall -Wextra -Wpedantic
NATIVE_LDLIBS = -lntl -lflint -lgmp -lzstd
NATIVE_BINARY := artifacts/native/census-native
NATIVE_SOURCES := $(wildcard native/*.cpp)
NATIVE_HEADERS := $(wildcard native/*.hpp)

CHECK_FILES := $(patsubst tests/%,check/%,$(filter-out tests/_%,$(wildcard tests/*.m))) \
               $(patsubst examples/%,check/%,$(filter-out examples/_%,$(wildcard examples/*.m)))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

census-native: $(NATIVE_BINARY)

$(NATIVE_BINARY): $(NATIVE_SOURCES) $(NATIVE_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(LDFLAGS) -o $@ $(NATIVE_SOURCES) $(NATIVE_LDLIBS) $(LDLIBS)

NATIVE_CHECKS := $(filter check/native_%,$(CHECK_FILES))
$(NATIVE_CHECKS) $(NATIVE_CHECKS:check/%=debug/%): census-native

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

.PHONY: all check check-no-dups census-native

