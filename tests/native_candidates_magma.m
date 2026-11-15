import "tests/_support.m": shell_quote;

// The standalone child returns its marker only after all assertions pass.
procedure test_native_candidates_magma()
    resultpath, result := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete result;
    delete logfile;
    try
        command := "\"${MAGMA:-magma}\" -bn resultpath:=" cat shell_quote(resultpath)
            cat " tests/native/candidates_magma.m > " cat shell_quote(logpath) cat " 2>&1";
        status := System(command);
        error if status ne 0 or Read(resultpath) ne "CANDIDATES_MAGMA_OK\n", Read(logpath);
    catch e
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_native_candidates_magma();
