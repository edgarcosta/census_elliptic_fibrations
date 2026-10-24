procedure expect_error(action, message)
    raised := false;
    try
        action();
    catch e
        assert message in Sprint(e);
        raised := true;
    end try;
    assert raised;
end procedure;

function shell_quote(s)
    return "'" cat Join(Split(s,"'"),"'\"'\"'") cat "'";
end function;
