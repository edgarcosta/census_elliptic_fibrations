import "Tests/_support.m": shell_quote;

// Every child uses the launcher's exit status and completion marker.
procedure reference_parallel_check(check, processes)
    assert System("mkdir -p artifacts/test-tmp") eq 0;
    scratch := GetCurrentDirectory() cat "/artifacts/test-tmp";
    files := [];
    command := "pids=''; ";
    for shard in [1..processes] do
        result,handle := TemporaryFile("reference-" : ParentDirectory := scratch); delete handle;
        log,handle := TemporaryFile("reference-" : ParentDirectory := scratch); delete handle;
        Append(~files,<result,log>);
        command cat:= Sprintf("REFERENCE_SHARD=%o REFERENCE_PROCESSES=%o ",
            shard,processes) cat Join([ shell_quote(s) : s in [
            "sh","scripts/run_reference_check.sh",check,result,log ] ]," ")
            cat " & pids=\"$pids $!\"; ";
    end for;
    command cat:= "status=0; for pid in $pids; do wait \"$pid\" || status=1; done; exit $status";
    status := System(command);
    for pair in files do
        if status ne 0 then printf "%o",Read(pair[2]); end if;
        Remove(pair[1]); Remove(pair[2]);
    end for;
    error if status ne 0, "MISMATCH: reference child check failed";
end procedure;
