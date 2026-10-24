// Run from the repository root. Select tests with target:= or filename:=.
// Assign exitsignal:= to return the number of failures; debug:= opens the debugger.
SetQuitOnError(true);
if assigned debug then
    SetDebugOnError(true);
end if;

// Load the current entry points here until the package spec is available.
load "ef_common.m";
load "ef3_models.m";
load "ef3_orbits.m";
load "k3_minimality_f3.m";
load "ef3_stream.m";
load "ef_orbits.m";
load "ef3_gmodule.m";
load "elliptic_fibration_orbits_f2.m";
load "k3_minimality_f2.m";
load "ef2_gmodule.m";
load "ef3_burnside.m";

if assigned filename then
    tests := [filename];
else
    tests := Split(Pipe("find Tests examples -maxdepth 1 -type f -name '*.m' | sort", ""), "\n");
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
end if;
if assigned exitsignal then
    exit #failed;
end if;
