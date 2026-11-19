SetQuitOnError(true);
SetNthreads(1);
SetColumns(0);
SetSeed(20261119);
inputs := eval Read(inputs);

function local_bound(symbol)
    s := Sprint(symbol);
    if s in ["II", "II*", "I0", "I1"] then return 0, 1; end if;
    if s eq "III" then return 1/2, 2; end if;
    if s eq "III*" then return 3/2, 2; end if;
    if s eq "IV" then return 2/3, 3; end if;
    if s eq "IV*" then return 4/3, 3; end if;
    if s[#s] eq "*" then
        n := StringToInteger(s[2..#s-1]);
        return 1 + n/4, 4;
    end if;
    n := StringToInteger(s[2..#s]);
    m := n div 2;
    return m*(n-m)/n, n;
end function;

for item in inputs do
    key, b, codes, points, gram, rank, torsion, torsion_points := Explode(item);
    R := PolynomialRing(GF(2));
    K<t> := FunctionField(GF(2));
    a := [K | R!Intseq(c, 2) : c in codes];
    E := EllipticCurve(a);
    basis := [E | [K!(R!Intseq(p[1], 2))/(R!Intseq(p[2], 2)),
                   K!(R!Intseq(p[3], 2))/(R!Intseq(p[4], 2))] : p in points];
    assert #basis eq rank and rank gt 0;
    G := HeightPairingMatrix(basis);
    assert G eq Matrix(Rationals(), rank, rank, &cat gram);
    T, mT := TorsionSubgroup(E);
    assert Invariants(T) eq torsion;
    given_torsion := [E | ];
    for p in torsion_points do
        if #p eq 0 then
            Append(~given_torsion, E!0);
        else
            Append(~given_torsion, E![K!(R!Intseq(p[1], 2))/(R!Intseq(p[2], 2)),
                                      K!(R!Intseq(p[3], 2))/(R!Intseq(p[4], 2))]);
        end if;
    end for;
    assert #given_torsion eq #T;
    assert Seqset(given_torsion) eq {mT(tors) : tors in T};
    Cmax := Rationals()!0;
    D := 1;
    for fibre in LocalInformation(E) do
        correction, denominator := local_bound(fibre[5]);
        Cmax +:= Degree(fibre[1])*correction;
        D := LCM(D, denominator);
    end for;
    scale := LCM([Denominator(x) : x in Eltseq(G)]);
    B := ChangeRing(LLLGram(ChangeRing(scale*G, Integers())), Rationals())/scale;
    norms := [Rationals() | ];
    for j in [1..rank] do
        assert B[j, j] gt 0;
        Append(~norms, B[j, j]);
        for i in [j+1..rank] do
            factor := B[i, j]/B[j, j];
            for k in [j+1..rank] do
                B[i, k] -:= factor*B[j, k];
            end for;
        end for;
    end for;
    mu := &+norms/4;
    assert mu lt 6 + 2*b - Cmax;
    N := D^rank*Determinant(G);
    assert IsIntegral(N);
    checks := 0;
    // Only primes with square divisibility can divide a lattice index.
    for ell in [f[1] : f in Factorization(Integers()!N) | f[2] ge 2] do
        for v in VectorSpace(GF(ell), rank) do
            if IsZero(v) then continue; end if;
            first := Min([j : j in [1..rank] | v[j] ne 0]);
            if v[first] ne 1 then continue; end if;
            P := &+[E | (Integers()!v[j])*basis[j] : j in [1..rank]];
            for tors in T do
                assert #DivisionPoints(P + mT(tors), ell) eq 0;
                checks +:= 1;
            end for;
        end for;
    end for;
    assert #DivisionPoints(2*basis[1], 2) gt 0;
    printf "AUDIT|%o|%o|%o|%o|%o\n", key, checks, mu, Cmax, D;
end for;
quit;
