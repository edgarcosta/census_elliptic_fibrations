import "tests/_support.m": shell_quote;

procedure test_native_tables_magma()
    resultpath, result := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete result;
    delete logfile;
    try
        command := "\"${MAGMA:-magma}\" -bn resultpath:=" cat shell_quote(resultpath)
            cat " tests/native/tables_magma.m > " cat shell_quote(logpath) cat " 2>&1";
        status := System(command);
        error if status ne 0 or Read(resultpath) ne "TABLES_MAGMA_OK\n", Read(logpath);
        printf "\n%o", Read(logpath);
    catch e
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_native_tables_magma();
