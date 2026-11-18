import "tests/support/_support.m": shell_quote;

// Compare the native command with package computations in an isolated child.
procedure test_native_model()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        Write(referencepath, "stale reference");
        command := "make census-native > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/support/magma/model_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/compare/model.py " cat shell_quote(referencepath)
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

test_native_model();
