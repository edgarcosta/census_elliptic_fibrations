import "Tests/_support.m": shell_quote;

// Pins reference checks without artifacts/, including failed or incomplete results.
procedure test_reference_result_protocol()
    directory := TemporaryDirectory("reference space ");
    launcher := directory cat "/scripts/run_reference_check.sh";
    scriptpath,script := TemporaryFile();
    resultpath,result := TemporaryFile();
    logpath,logfile := TemporaryFile();
    delete script;
    Remove(scriptpath);
    scriptpath cat:= " input.m";
    delete result;
    delete logfile;
    try
        assert System("mkdir " cat shell_quote(directory cat "/scripts")) eq 0;
        assert System("cp scripts/run_reference_check.sh " cat shell_quote(launcher)) eq 0;
        assert System("test ! -e " cat shell_quote(directory cat "/artifacts")) eq 0;
        for body in [ "Write(resultpath,\"REFERENCE_CHECK_OK\");\nexit 0;\n",
            "assert false;\nWrite(resultpath,\"REFERENCE_CHECK_OK\");\n",
            "this is invalid magma;\n", "exit 0;\n" ] do
            Write(scriptpath,body : Overwrite := true);
            empty := Open(resultpath,"w");
            delete empty;
            command := Join([ shell_quote(s) : s in [ "sh",launcher,scriptpath,resultpath,logpath ] ]," ");
            status := System(command);
            if body[1] eq "W" then
                error if status ne 0,Read(logpath);
                assert Read(resultpath) eq "REFERENCE_CHECK_OK\n";
            else
                assert status ne 0 and Read(resultpath) ne "REFERENCE_CHECK_OK\n";
            end if;
        end for;
    catch e
        Remove(scriptpath);
        Remove(resultpath);
        Remove(logpath);
        assert System("rm -rf -- " cat shell_quote(directory)) eq 0;
        error e;
    end try;
    Remove(scriptpath);
    Remove(resultpath);
    Remove(logpath);
    assert System("rm -rf -- " cat shell_quote(directory)) eq 0;
end procedure;

test_reference_result_protocol();
