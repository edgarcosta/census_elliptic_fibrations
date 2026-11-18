import "tests/support/_support.m": shell_quote;

// The native CLI must build and reject invalid invocations without arithmetic.
procedure test_native_build()
    script := GetCurrentDirectory() cat "/tests/native/build.py";
    assert System("python3 " cat shell_quote(script)) eq 0;
end procedure;

test_native_build();
