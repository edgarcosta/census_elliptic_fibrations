SetQuitOnError(true);
SetNthreads(1);
SetColumns(0);
AttachSpec("spec");
if not assigned prefixes then prefixes := ""; end if;
if not assigned prefix then prefix := ""; end if;

q := StringToInteger(field);
R<s, t> := PolynomialRing(GF(q), 2);
if action eq "count" then
    indices := [StringToInteger(x) : x in Split(prefixes, ",")];
    P := FibrationPrefixes(R);
    n, cells, rows := NumberOfFibrations(R : PrefixIndices := indices,
        GenericFibreKinds := ["smooth"], Minimalities := ["everywhere minimal"]);
    file := Open(destination, "w");
    fprintf file, "prefix,count,form_1_degree,form_1_key,form_2_degree,form_2_key\n";
    degrees := q eq 2 select [2, 6] else [4, 8];
    for row in rows do
        i := row[1];
        keys := [Seqint([Integers() | MonomialCoefficient(P[i][j], s^k*t^(degrees[j]-k))
            : k in [0..degrees[j]]], q) : j in [1, 2]];
        fprintf file, "%o,%o,%o,%o,%o,%o\n", i,
            row[2][<"smooth", "everywhere minimal">], degrees[1], keys[1], degrees[2], keys[2];
    end for;
    delete file;
elif action eq "write" then
    n := WriteFibrationRepresentativeKeys(R, destination :
        PrefixIndices := [StringToInteger(prefix)], GenericFibreKinds := ["smooth"],
        Minimalities := ["everywhere minimal"]);
else
    error "Unknown key driver action";
end if;
quit;
