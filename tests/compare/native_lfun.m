import "tests/support/_support.m": compare_native;
import "tests/support/magma/lfun_native_reference.m": write_lfun_reference;

compare_native(write_lfun_reference, "tests/compare/lfun.py");
