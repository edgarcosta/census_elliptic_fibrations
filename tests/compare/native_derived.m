import "tests/support/_support.m": compare_native;
import "tests/support/magma/derived_reference.m": write_reference;

compare_native(write_reference, "tests/compare/derived.py");
