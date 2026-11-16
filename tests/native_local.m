import "tests/_support.m": shell_quote;

// Compare every local field and sign through the native public command.
procedure test_native_local()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        Write(referencepath, "stale reference");
        command := "make census-native > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/native/local_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/native/local.py " cat shell_quote(referencepath)
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

test_native_local();
