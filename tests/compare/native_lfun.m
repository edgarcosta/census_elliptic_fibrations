import "tests/support/_support.m": shell_quote;

// Completed signed coefficients must agree with the package in both charts.
procedure test_native_lfun()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        Write(referencepath, "stale reference");
        command := "make census-native census-native-scalar > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/support/magma/lfun_native_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/compare/lfun.py " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        Remove(referencepath);
        Remove(logpath);
        error e;
    end try;
    Remove(referencepath);
    Remove(logpath);
end procedure;

test_native_lfun();
