import "tests/support/_support.m": shell_quote;

procedure test_accessible_f2_bases()
    script := GetCurrentDirectory() cat "/tests/support/search_f2.py";
    assert System("python3 " cat shell_quote(script)) eq 0;
end procedure;

test_accessible_f2_bases();
