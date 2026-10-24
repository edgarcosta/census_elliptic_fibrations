// Test and timing for one complete dimension-13 fiber, without disk output.
load "ef_common.m";
load "ef3_models.m";
load "ef3_orbits.m";
load "k3_minimality_f3.m";
load "ef3_stream.m";

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
quit;
