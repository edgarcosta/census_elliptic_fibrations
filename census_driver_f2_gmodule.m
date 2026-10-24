// Prints, for each of the 88 prefixes <a1,a3>, the fibre group order, the
// fibre dimension and the orbit counts, then the census totals and CPU time.
// Run magma -b census_driver_f2_gmodule.m from the repository root.

load "ef_common.m";
load "ef_orbits.m";
load "ef2_gmodule.m";

R<s,t> := PolynomialRing(GF(2),2);
start := Cputime();
prefixes := FibrationOrbitPrefixesF2G(R);
all_orbits := 0; nonsingular := 0; k3 := 0;
print "i, a1, a3, #E, dim, all, Delta != 0, K3, CPU s";
for i in [1..#prefixes] do
    clock := Cputime();
    P := prefixes[i];
    F := FibrationFibreF2G(P);
    c := EF_FibreCounts(F);
    all_orbits +:= c[1]; nonsingular +:= c[2]; k3 +:= c[3];
    printf "%o, %o, %o, %o, %o, %o, %o, %o, %o\n", i, P[1], P[2], #F`Group,
        F`Dimension, c[1], c[2], c[3], Cputime(clock);
end for;
assert all_orbits eq 1402468 and nonsingular eq 1401732;
printf "All orbits: %o\nDelta != 0 orbits: %o\nK3 orbits: %o\n",
    all_orbits, nonsingular, k3;
printf "Total CPU seconds: %o\n", Cputime(start);
quit;
