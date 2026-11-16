import "tests/_support.m": shell_quote;

// Fresh package tables and direct counts check the native counting command.
procedure test_native_tables()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        Write(referencepath, "stale reference");
        command := "make census-native > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/native/tables_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/native/tables.py " cat shell_quote(referencepath)
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

test_native_tables();
