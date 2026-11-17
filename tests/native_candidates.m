import "tests/_support.m": shell_quote;

procedure test_native_candidates()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    resultpath, result := TemporaryFile();
    delete reference;
    delete logfile;
    delete result;
    try
        command := "make census-native artifacts/native/candidates-exact > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/native/candidates_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath cat ".supplied")
            cat " resultpath:=" cat shell_quote(resultpath)
            cat " tests/native/candidates_magma.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0 or Read(resultpath) ne "CANDIDATES_MAGMA_OK\n", Read(logpath);
        command := "python3 tests/native/candidates.py " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        Remove(referencepath);
        if FileExists(referencepath cat ".supplied") then Remove(referencepath cat ".supplied"); end if;
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(referencepath);
    Remove(referencepath cat ".supplied");
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_native_candidates();
