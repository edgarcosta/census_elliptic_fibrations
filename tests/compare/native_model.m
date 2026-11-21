import "tests/support/_support.m": compare_native;
import "tests/support/magma/model_reference.m": write_model_reference;

compare_native(write_model_reference, "tests/compare/model.py");
