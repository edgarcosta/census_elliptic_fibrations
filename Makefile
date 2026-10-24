MAGMA ?= magma

CHECK_FILES := $(patsubst Tests/%,check/%,$(wildcard Tests/*.m)) \
               $(patsubst examples/%,check/%,$(wildcard examples/*.m))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

check/%.m: Tests/%.m | check-no-dups
	@$(MAGMA) -b filename:=$< exitsignal:='' run_tests.m
check/%.m: examples/%.m | check-no-dups
	@$(MAGMA) -b filename:=$< exitsignal:='' run_tests.m

debug/%.m: Tests/%.m | check-no-dups
	@$(MAGMA) -b filename:=$< debug:='' run_tests.m
debug/%.m: examples/%.m | check-no-dups
	@$(MAGMA) -b filename:=$< debug:='' run_tests.m

check: $(CHECK_FILES)

check-no-dups:
	@sh scripts/no_dup_test_basenames.sh

.PHONY: all check check-no-dups
