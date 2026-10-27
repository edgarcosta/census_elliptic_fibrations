import "Tests/_support.m": shell_quote;

// Pins absolute attachment in a fresh process whose working directory is elsewhere.
procedure test_fresh_attachment()
    root := GetCurrentDirectory();
    directory := TemporaryDirectory("attachment space ");
    resultpath,result := TemporaryFile();
    logpath,logfile := TemporaryFile();
    delete result;
    delete logfile;
    try
        status := System(Join([ shell_quote(s) : s in [ "sh",root cat "/scripts/magma_child.sh",
            root cat "/Tests/_attachment_child.m",directory,root cat "/spec",resultpath,logpath ] ]," "));
        error if status ne 0 or Read(resultpath) ne "ATTACHMENT_OK\n",Read(logpath);
    catch e
        Remove(resultpath);
        Remove(logpath);
        assert System("rmdir " cat shell_quote(directory)) eq 0;
        error e;
    end try;
    Remove(resultpath);
    Remove(logpath);
    assert System("rmdir " cat shell_quote(directory)) eq 0;
end procedure;

test_fresh_attachment();
