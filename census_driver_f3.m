// Prints, for each of the 1,860 prefixes <C,A>, the stabilizer order, the
// fibre dimension and the orbit counts, then the census totals and CPU time.
// Run magma -b census_driver_f3.m from the repository root.

load "ef_common.m";
load "ef3_models.m";
load "ef_orbits.m";
load "ef3_gmodule.m";

R<s,t> := PolynomialRing(GF(3),2);
start := Cputime();
prefixes := FibrationOrbitPrefixesF3G(R);
all_orbits := 0; nonsingular := 0; k3 := 0; k3_c0 := 0;
print "i, C, A, #H, dim, all, Delta != 0, K3, CPU s";
for i in [1..#prefixes] do
    clock := Cputime();
    P := prefixes[i];
    F := FibrationFibreF3G(P);
    d := F`Dimension;
    c := EF_FibreCounts(F);
    all_orbits +:= c[1]; nonsingular +:= c[2]; k3 +:= c[3];
    if P[1] eq 0 then k3_c0 +:= c[3]; end if;
    printf "%o, %o, %o, %o, %o, %o, %o, %o, %o\n", i, P[1], P[2], #P[3], d,
        c[1], c[2], c[3], Cputime(clock);
end for;
assert all_orbits eq 1307814237 and nonsingular eq 1307813899;
printf "All orbits: %o\nDelta != 0 orbits: %o\nK3 orbits: %o\n",
    all_orbits, nonsingular, k3;
printf "K3 with C = 0: %o\nK3 with C != 0: %o\n", k3_c0, k3-k3_c0;
printf "Non-minimal somewhere: %o\nTotal CPU seconds: %o\n",
    nonsingular-k3, Cputime(start);
quit;
