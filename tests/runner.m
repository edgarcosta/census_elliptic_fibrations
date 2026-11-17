import "tests/_support.m": expect_error, shell_quote;

procedure test_runner_exit_status()
    magma := GetEnv("MAGMA");
    if magma eq "" then magma := "magma"; end if;
    fragment, handle := TemporaryFile();
    delete handle;
    logpath, handle := TemporaryFile();
    delete handle;
    try
        for example in [<"assert true;", true>, <"assert false;", false>] do
            Remove(fragment);
            Write(fragment, example[1]);
            command := Join([shell_quote(s) : s in [magma, "-bn",
                "filename:=" cat fragment, "run_tests.m"]], " ");
            status := System(command cat " > " cat shell_quote(logpath) cat " 2>&1");
            error if (status eq 0) ne example[2], Sprintf(
                "Runner exit status %o for %o:\n%o", status, example[1], Read(logpath));
        end for;
    catch e
        Remove(fragment);
        Remove(logpath);
        error e;
    end try;
    Remove(fragment);
    Remove(logpath);
end procedure;

test_runner_exit_status();

// A traceback can contain the expected text as an argument to expect_error.
procedure test_expected_error_message()
    raised := false;
    try
        expect_error(procedure() x := 1 div 0; end procedure,
            "expected error");
    catch e
        raised := true;
    end try;
    assert raised;
end procedure;

test_expected_error_message();
