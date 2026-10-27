// Times one complete prefix of each kind, then one of the commonest kind with
// representatives written. Run from the repository root; it writes
// bench_f3_nowrite* and bench_f3_write* files in the current directory.
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";

R<s,t> := PolynomialRing(GF(3),2);
P := FibrationOrbitPrefixesF3(R);
buckets := FibrationBenchmarkBucketsF3(P);
keys := Sort([b : b in Keys(buckets)]);
for b in keys do
    printf "bucket %o: %o prefixes\n", b, #buckets[b];
end for;
indices := [Min(buckets[b]) : b in keys];
printf "bucket representatives: %o\n", indices;
results := CensusF3(R,indices,"bench_f3_nowrite" : WriteRepresentatives:=false);
for b in keys do
    i := Min(buckets[b]);
    st := [z[2] : z in results | z[1] eq i][1];
    printf "TIME %o %o %o\n", b, #buckets[b], st`CPUSeconds;
end for;
dominant := [j : j in buckets[<false,13,1,true>] | P[j][2] ne 0][1];
w := CensusF3(R,[dominant],"bench_f3_write");
printf "WRITE %o %o\n", dominant, w[1][2]`CPUSeconds;
quit;
