import "tests/support/_support.m": shell_quote;

// Pins recovery and receipt rejection through the production commands.
procedure test_stage1_commands()
    logpath := "artifacts/stage1-tests.log";
    command := "make -j2 census-native > " cat shell_quote(logpath) cat " 2>&1";
    error if System(command) ne 0, Read(logpath);
    command := "python3 tests/compare/run_stage1.py >> " cat shell_quote(logpath) cat " 2>&1";
    error if System(command) ne 0, Read(logpath);
end procedure;

test_stage1_commands();
