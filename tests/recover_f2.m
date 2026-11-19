import "tests/support/_support.m": shell_quote;

// Retained deeper bases need a checked import; failed or misassociated work stays open.
procedure test_recovered_f2_bases()
    script := GetCurrentDirectory() cat "/tests/support/recover_f2.py";
    assert System("python3 " cat shell_quote(script)) eq 0;
end procedure;

test_recovered_f2_bases();
