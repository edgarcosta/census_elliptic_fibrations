SetQuitOnError(true);
// Run from the repository root. Select tests with target:= or filename:=.
// Assign debug:= to open the debugger on failure.
if assigned debug then
    SetQuitOnError(false);
    SetDebugOnError(true);
end if;

load "config.m";

if assigned filename then
    tests := [filename];
else
    tests := Split(Pipe("python3 tests/list_tests.py", ""), "\n");
end if;
if not assigned target then
    target := "";
end if;

failed := [];
counter := 0;
for filename in tests do
    if target in filename then
        counter +:= 1;
        timestamp := Time();
        printf "%o: ", filename;
        if assigned debug then
            assert eval (Read(filename) cat "return true;");
            printf "Success! %o s\n", Time(timestamp);
        else
            try
                assert eval (Read(filename) cat "return true;");
                printf "Success! %o s\n", Time(timestamp);
            catch e
                Append(~failed, filename);
                printf "Fail! %o s\n%o\n", Time(timestamp), e;
            end try;
        end if;
    end if;
end for;
if counter eq 0 then
    print "No matching target";
    exit 1;
end if;
if #failed gt 0 then
    print "Tests failed:";
    for filename in failed do
        print filename;
    end for;
    exit #failed;
end if;
quit;
