import "tests/support/_archive.m": test_archive;
import "tests/support/_archive_oracle.m": test_archive_canonicalization;

// All F_2 prefixes compare complete orbit multisets and archived classification.
procedure test_archive_fixtures_f2()
    test_archive_canonicalization(2);
    test_archive(2);
end procedure;

test_archive_fixtures_f2();
