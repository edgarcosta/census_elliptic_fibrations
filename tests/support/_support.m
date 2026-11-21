procedure expect_error(action, message)
    raised := false;
    try
        action();
    catch e
        assert message in Sprint(e`Object);
        raised := true;
    end try;
    assert raised;
end procedure;

function shell_quote(s)
    return "'" cat Join(Split(s, "'"), "'\"'\"'") cat "'";
end function;

procedure compare_native(write_reference, script)
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        command := "make census-native > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        write_reference(referencepath);
        command := "python3 " cat shell_quote(script) cat " " cat shell_quote(referencepath)
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
