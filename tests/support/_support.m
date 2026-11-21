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

procedure compare_native(write_reference, script : Targets := "", Suffixes := [])
    referencepath, reference := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete reference;
    delete logfile;
    try
        command := "make census-native " cat Targets cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
        write_reference(referencepath);
        command := "python3 " cat shell_quote(script) cat " " cat shell_quote(referencepath)
            cat " > " cat shell_quote(logpath) cat " 2>&1";
        error if System(command) ne 0, Read(logpath);
    catch e
        for suffix in [""] cat Suffixes do
            if FileExists(referencepath cat suffix) then Remove(referencepath cat suffix); end if;
        end for;
        Remove(logpath);
        error e;
    end try;
    for suffix in [""] cat Suffixes do
        if FileExists(referencepath cat suffix) then Remove(referencepath cat suffix); end if;
    end for;
    Remove(logpath);
end procedure;
