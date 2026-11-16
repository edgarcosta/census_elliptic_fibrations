import "tests/_support.m": shell_quote;

procedure test_native_records()
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    resultpath, result := TemporaryFile();
    delete reference;
    delete logfile;
    delete result;
    try
        command := "make census-native artifacts/native/records-exact > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn referencepath:=" cat shell_quote(referencepath)
            cat " tests/native/records_reference.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        command := "\"${MAGMA:-magma}\" -bn codecpath:=" cat shell_quote(referencepath cat ".codec")
            cat " resultpath:=" cat shell_quote(resultpath)
            cat " tests/native/records_magma.m > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0 or Read(resultpath) ne "RECORDS_MAGMA_OK\n", Read(logpath);
        command := "python3 tests/native/records.py " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        Remove(referencepath);
        Remove(referencepath cat ".codec");
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(referencepath);
    Remove(referencepath cat ".codec");
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_native_records();
