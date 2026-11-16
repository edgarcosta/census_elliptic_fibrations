import "tests/_support.m": shell_quote;

// Exact fields agree on both charts, including the supplied constant-field controls.
procedure test_native_derived()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        Write(referencepath, "stale reference");
        command := "make census-native artifacts/native/derived-exact > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/native/derived_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath cat ".exact")
            cat " tests/native/derived_exact.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "python3 tests/native/derived.py " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        Remove(referencepath);
        Remove(referencepath cat ".exact");
        Remove(logpath);
        error e;
    end try;
    Remove(referencepath);
    Remove(referencepath cat ".exact");
    Remove(logpath);
end procedure;

test_native_derived();
