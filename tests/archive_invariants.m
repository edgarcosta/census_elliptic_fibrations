import "tests/support/_archive.m": test_archive;
import "tests/support/_archive_oracle.m": test_archive_canonicalization;

// Two complete F_3 prefixes exercise both canonicalization branches.
procedure test_archive_invariants_f3()
    test_archive_canonicalization(3);
    test_archive(3);
end procedure;

test_archive_invariants_f3();
