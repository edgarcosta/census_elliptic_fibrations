import "tests/support/_weierstrass_change.m": equation_change;

function sample_form(R, mask, degree)
    return &+[R | ((mask div 2^j) mod 2)*R.1^j*R.2^(degree-j)
        : j in [0..degree]];
end function;

// Between normalized F_2 models, a2 = 0 forces r = q^2+a1*q.
// Search only this row's coordinate orbit against public prefix listings.
procedure match_sample_model(U, prefixes, ~listings, ~key)
    R<s, t> := Universe(U);
    for mask in [0..7] do
        q := sample_form(R, mask, 2);
        B := equation_change(U, GF(2)!1, q^2+U[1]*q, q, R!0);
        for g in GL(2, 2) do
            sub := [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t];
            A := [Evaluate(a, sub) : a in B];
            assert A[2] eq 0;
            i := Index(prefixes, <A[1], A[3]>);
            if i eq 0 then continue; end if;
            if not IsDefined(listings, i) then
                listings[i] := Seqset(FibrationRepresentativeKeys(R, i :
                    GenericFibreKinds := ["smooth"], Minimalities := ["everywhere minimal"]));
            end if;
            for vmask in [0..127] do
                v := sample_form(R, vmask, 6);
                V := [R | A[1], 0, A[3], A[4]+v*A[1], A[5]+v*A[3]+v^2];
                candidate := FibrationKey(V);
                if candidate in listings[i] then
                    assert equation_change(A, GF(2)!1, R!0, R!0, v) eq V;
                    assert aInvariantsFromKey(R, candidate) eq V;
                    key := candidate;
                    return;
                end if;
            end for;
        end for;
    end for;
    error "No listed representative in the coordinate orbit of", FibrationKey(U);
end procedure;

// These invariants do not depend on the chosen orbit representative or bad-place names.
procedure check_sample_invariants(W, row)
    F<T> := FunctionField(GF(2));
    E := EllipticCurve(W, F : Infinity := true);
    G := TorsionSubgroup(E);
    torsion := IsEmpty(Invariants(G)) select "trivial"
        else Join(["Z/" cat Sprint(n) : n in Invariants(G)], " x ");
    assert torsion eq row[4];
    fibres := {**};
    for info in LocalInformation(E) do
        Include(~fibres, Sprintf("%o:%o:%o:%o", Degree(info[1]), info[3], info[4], info[5]));
    end for;
    assert fibres eq SequenceToMultiset(Split(row[5], ";"));
    L := LFunction(E);
    assert Coefficients(L) eq [StringToInteger(x) : x in Split(row[6], " ")];
    assert LeadingCoefficient(L) eq StringToInteger(row[7])*2^Degree(L);
    analytic := AnalyticInformation(E, L);
    assert analytic[1] eq StringToInteger(row[8]);
end procedure;

// The published F_2 census table is incomplete; this sample cannot check whole-census cell counts.
procedure test_published_f2()
    R<s, t> := PolynomialRing(GF(2), 2);
    lines := Split(Read("tests/fixtures/published-f2/sample.csv"), "\n");
    assert lines[1] eq "source_row,key,a1,mw_torsion,singular_fibres,L_function,L_sign,analytic_rank";
    assert #lines eq 17;
    prefixes := FibrationPrefixes(R);
    listings := AssociativeArray();
    classes := AssociativeArray();
    changed := false;
    previous := 0;
    for line in lines[2..#lines] do
        row := Split(line, ",");
        source_row := StringToInteger(row[1]);
        assert source_row gt previous;
        previous := source_row;
        table_key := StringToInteger(row[2]);
        U := aInvariantsFromKey(R, table_key);
        assert Sprint(U[1]) eq row[3];
        if not IsDefined(classes, row[3]) then classes[row[3]] := 0; end if;
        classes[row[3]] +:= 1;
        match_sample_model(U, prefixes, ~listings, ~key);
        changed or:= key ne table_key;
        V := aInvariantsFromKey(R, key);
        for W in [U, V] do
            kind, minimality := FibrationProperties(W);
            assert <kind, minimality> eq <"smooth", "everywhere minimal">;
            check_sample_invariants(W, row);
        end for;
    end for;
    assert #Keys(classes) eq 4 and &and[n eq 4 : n in Values(classes)];
    assert changed;
end procedure;

test_published_f2();
