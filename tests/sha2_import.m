procedure test_sha2_evidence_import()
    status := System("mkdir -p artifacts && python3 tests/support/sha2_import.py > artifacts/sha2-import-tests.log 2>&1");
    error if status ne 0, Read("artifacts/sha2-import-tests.log");
end procedure;

test_sha2_evidence_import();
