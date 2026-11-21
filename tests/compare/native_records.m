import "tests/support/_support.m": compare_native;
import "tests/support/magma/records_reference.m": write_records_reference;

compare_native(write_records_reference, "tests/compare/records.py");
