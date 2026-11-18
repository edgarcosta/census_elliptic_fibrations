import "tests/support/_support.m": shell_quote;

procedure test_native_candidates()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        command := "make census-native artifacts/native/candidates-exact > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/support/magma/candidates_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/compare/candidates.py " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        Remove(referencepath);
        if FileExists(referencepath cat ".supplied") then Remove(referencepath cat ".supplied"); end if;
        Remove(logpath);
        error e;
    end try;
    Remove(referencepath);
    Remove(referencepath cat ".supplied");
    Remove(logpath);
end procedure;

test_native_candidates();
