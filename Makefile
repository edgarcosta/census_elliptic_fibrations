MAGMA ?= magma
MAGMA_PATH := $(shell command -v "$(MAGMA)")

CXXFLAGS ?= -O2 -Wall -Wextra -Wpedantic
NATIVE_LDLIBS = -lntl -lflint -lgmpxx -lgmp
NATIVE_BINARY := artifacts/native/census-native

TEST_PATHS := $(shell python3 tests/list_tests.py)
ifneq ($(.SHELLSTATUS),0)
$(error tests/list_tests.py failed)
endif
CHECK_FILES := $(addprefix check/,$(notdir $(TEST_PATHS)))

MAKEFLAGS += -j$(shell nproc 2>/dev/null || printf 1)

all: check

census-native: $(NATIVE_BINARY)

NATIVE_SOURCES := native/main.cpp $(wildcard native/*/*.cpp)
NATIVE_HEADERS := $(wildcard native/*/*.hpp)

$(NATIVE_BINARY): $(NATIVE_SOURCES) $(NATIVE_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(LDFLAGS) -o $@ $(NATIVE_SOURCES) $(NATIVE_LDLIBS) $(LDLIBS)

NATIVE_CHECKS := $(filter check/native_%,$(CHECK_FILES))
$(NATIVE_CHECKS) $(NATIVE_CHECKS:check/%=debug/%): census-native

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

.PHONY: all check census-native

artifacts/native/candidates-exact: tests/native/candidates_exact.cpp native/invariants/invariants.cpp $(NATIVE_HEADERS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(LDFLAGS) -o $@ tests/native/candidates_exact.cpp native/invariants/invariants.cpp $(NATIVE_LDLIBS) $(LDLIBS)
