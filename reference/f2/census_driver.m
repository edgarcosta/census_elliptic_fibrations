// Prints, for each of the 88 prefixes <a1,a3>, the orbit counts and the number
// of everywhere-minimal (K3) models, then the census totals and CPU times.
// Run magma -b reference/f2/census_driver.m from the repository root.

load "reference/f2/elliptic_fibration_orbits_f2.m";
load "reference/f2/k3_minimality_f2.m";

R<s,t> := PolynomialRing(GF(2),2);
start := Cputime();
prefixes := FibrationOrbitPrefixesF2(R);
a1_reps := [R|0,s^2,s*t,s^2+s*t+t^2];
k3_by_a1 := [0,0,0,0];
all_orbits := 0; nonsingular := 0; k3 := 0;
one_place := 0; several_places := 0;
enumeration_time := 0.0; filter_time := 0.0;
print "i, a1, <dim,#maps,#orbits,#retained>, K3, enumeration CPU s, filter CPU s";
for i in [1..#prefixes] do
    clock := Cputime();
    K,I := FibrationFiberOrbitKeysF2(prefixes[i]);
    enumeration_seconds := Cputime(clock);
    clock := Cputime();
    count := 0;
    for key in K do
        A := FibrationFromKeyF2(R,key);
        ok,bad := IsEverywhereMinimalF2(A);
        if ok then
            count +:= 1;
        elif #bad eq 1 then
            one_place +:= 1;
        elif #bad ge 2 then
            several_places +:= 1;
        end if;
    end for;
    filter_seconds := Cputime(clock);
    assert #K eq I[4];
    a1 := prefixes[i][1];
    group := Index(a1_reps,a1);
    assert group ne 0;
    k3_by_a1[group] +:= count;
    all_orbits +:= I[3]; nonsingular +:= I[4]; k3 +:= count;
    enumeration_time +:= enumeration_seconds;
    filter_time +:= filter_seconds;
    printf "%o, %o, %o, %o, %o, %o\n", i,a1,I,count,
        enumeration_seconds,filter_seconds;
end for;
assert nonsingular eq k3+one_place+several_places;
printf "All orbits: %o\nDelta != 0 orbits: %o\nK3 orbits: %o\n",
    all_orbits,nonsingular,k3;
for j in [1..4] do
    printf "K3 with a1 = %o: %o\n", a1_reps[j],k3_by_a1[j];
end for;
printf "Non-minimal at exactly 1 place: %o\n", one_place;
printf "Non-minimal at 2 or more places: %o\n", several_places;
printf "Enumeration CPU seconds: %o\nFilter CPU seconds: %o\nTotal CPU seconds: %o\n",
    enumeration_time,filter_time,Cputime(start);
quit;
