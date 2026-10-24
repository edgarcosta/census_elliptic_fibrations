// Run from the repository root: magma -b filename:=examples/slow/check_burnside.m exitsignal:='' run_tests.m
// The constants below were computed independently in Python.

R<s,t> := PolynomialRing(GF(3),2);
clock := Cputime();
EFB_CheckActionConventions(R);
n,spectrum,empty,fixedSum,G := EFB_BurnsideF3(R);
prefixCounts := EFB_PrefixBurnsideF3(R);
printf "Prefix counts <quartic C, C!=0, C=0>: %o\n", prefixCounts;
printf "Effective group order: %o\n", Order(G);
print "Fixed-locus dimension, number of group elements";
for d in [0..27] do
    if spectrum[d] ne 0 then printf "%o, %o\n", d,spectrum[d]; end if;
end for;
printf "Empty fixed loci: %o\nSum of fixed points: %o\n", empty,fixedSum;
printf "All model orbits: %o\nCPU seconds: %o\n", n,Cputime(clock);

// Optional second enumeration of the same effective group:
// n2,sp2,e2,sum2,G2 := EFB_BurnsideF3(R : UseGroupIterator := true);
// assert n2 eq n and e2 eq empty and sum2 eq fixedSum;
// for d in [0..27] do assert sp2[d] eq spectrum[d]; end for;
