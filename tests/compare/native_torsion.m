import "tests/support/_support.m": compare_native;
import "tests/support/magma/torsion_reference.m": write_torsion_reference;

compare_native(write_torsion_reference, "tests/compare/torsion.py");
