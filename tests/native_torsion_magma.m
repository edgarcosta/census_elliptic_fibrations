import "tests/_support.m": shell_quote;

// The child has an explicit quit, so it runs outside the shared test runner.
procedure test_native_torsion_magma()
    resultpath, result := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete result;
    delete logfile;
    try
        command := "\"${MAGMA:-magma}\" -bn resultpath:=" cat shell_quote(resultpath)
            cat " tests/native/torsion_magma.m > " cat shell_quote(logpath) cat " 2>&1";
        status := System(command);
        error if status ne 0 or Read(resultpath) ne "TORSION_MAGMA_OK\n", Read(logpath);
    catch e
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_native_torsion_magma();
