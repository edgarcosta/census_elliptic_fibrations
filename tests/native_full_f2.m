import "tests/_support.m": shell_quote;

// Checks every F_2 orbit and shared published field, with fixed F_3 controls.
procedure test_native_full_f2()
    logpath := "artifacts/full-f2/check.log";
    assert System("mkdir -p artifacts/full-f2") eq 0;
    command := "make census-native artifacts/native/f2-orbit-match > "
        cat shell_quote(logpath) cat " 2>&1";
    error if System(command) ne 0, Read(logpath);
    command := "python3 tests/native/full_f2.py >> " cat shell_quote(logpath) cat " 2>&1";
    error if System(command) ne 0, Read(logpath);
end procedure;

test_native_full_f2();
