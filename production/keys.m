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
elif action eq "controls" then
    keys := q eq 2 select [4884480, 11821184] else [2729765910770, 66249649205, 320009441124];
    extra := q eq 2 select [[R | s*t, 0, s^6+t^6, 0, 0],
        [R | 0, 0, t^6, s^7*t, 0]] else [[R | 0, s^4+t^4, 0, s^4*t^4, 0]];
    keys cat:= [FibrationKey(W) : W in extra];
    Write(destination, FibrationInvariantCSV());
    F<u> := FunctionField(GF(q));
    Z<T> := PolynomialRing(Integers());
    for key in keys do
        W := aInvariantsFromKey(R, key);
        row := FibrationInvariantRow(key, W);
        assert row`Status eq "ok";
        E := EllipticCurve([F | Evaluate(a, [u, 1]) : a in W]);
        L := Z!([1] cat row`LowerL);
        for j in [0..((row`Degree-1) div 2)] do
            L +:= row`Sign*q^(row`Degree-2*j)*Coefficient(L, j)*T^(row`Degree-j);
        end for;
        assert L eq Z!LFunction(E);
        if key eq 2729765910770 then assert row`AdmissibleSha eq [1, 4, 9]; end if;
        if key eq 66249649205 then assert row`AdmissibleSha eq [1]; end if;
        Write(destination, FibrationInvariantCSV(row));
    end for;
else
    error "Unknown key driver action";
end if;
quit;
