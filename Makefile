MAGMA ?= magma
MAGMA_PATH := $(shell command -v "$(MAGMA)")

TEST_PATHS := $(shell python3 tests/list_tests.py)
ifneq ($(.SHELLSTATUS),0)
$(error tests/list_tests.py failed)
endif
CHECK_FILES := $(addprefix check/,$(notdir $(TEST_PATHS)))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

define test_rule
check/$(notdir $(1)): $(1)
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$(1) tests/run_tests.m
debug/$(notdir $(1)): $(1)
	@mkdir -p artifacts/test-tmp
	@MAGMA="$(MAGMA_PATH)" MAGMA_TEMP_DIR="$(CURDIR)/artifacts/test-tmp" "$(MAGMA)" -bn filename:=$(1) debug:='' tests/run_tests.m
endef
$(foreach file,$(TEST_PATHS),$(eval $(call test_rule,$(file))))

check: $(CHECK_FILES)

check-no-dups:
	@python3 tests/list_tests.py > /dev/null

.PHONY: all check check-no-dups
