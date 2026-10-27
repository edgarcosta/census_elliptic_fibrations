import "Tests/_support.m": shell_quote;

// Runs the independent check in a fresh process and verifies its result file.
procedure test_reference_burnside()
    resultpath, result := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete result;
    delete logfile;
    try
        command := Join([ shell_quote(s) : s in [ "sh",
            "scripts/run_reference_check.sh", "reference/checks/check_burnside.m",
            resultpath, logpath ] ], " ");
        status := System(command);
        error if status ne 0, Read(logpath);
        error if Read(resultpath) ne "REFERENCE_CHECK_OK\n", Read(logpath);
    catch e
        Remove(resultpath);
        Remove(logpath);
        error e;
    end try;
    Remove(resultpath);
    Remove(logpath);
end procedure;

test_reference_burnside();
