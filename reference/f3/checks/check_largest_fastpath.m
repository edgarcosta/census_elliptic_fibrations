// Test and timing for one complete dimension-13 fiber, without disk output.
load "reference/f3/ef_common.m";
load "reference/f3/ef3_models.m";
load "reference/f3/ef3_orbits.m";
load "reference/f3/k3_minimality_f3.m";
load "reference/f3/ef3_stream.m";

procedure EF3_CheckLargestFastpath()
    R<s,t> := PolynomialRing(GF(3),2);
    P := FibrationOrbitPrefixesF3(R);
    buckets := FibrationBenchmarkBucketsF3(P);
    i := [j : j in buckets[<false,13,1,true>] | P[j][2] ne 0][1];
    start := Cputime();
    stats := FibrationStreamPrefixF3(P[i],"unused" : WriteRepresentatives:=false);
    assert stats`Complete and stats`K3Orbits eq 3^13;
    printf "Largest trivial automatic prefix %o: %o K3; %o CPU s; no file I/O\n",
        i,stats`K3Orbits,Cputime(start);
end procedure;

EF3_CheckLargestFastpath();
if assigned resultpath then
    Write(resultpath,"REFERENCE_CHECK_OK");
end if;
quit 0;
