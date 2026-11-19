import "_archive_invariant.m": invariant_context, orbit_invariant;
import "_weierstrass_change.m": equation_change;

// Direct equation substitutions check the optimized canonical form independently.
procedure test_archive_canonicalization(q)
    R<s, t> := PolynomialRing(GF(q), 2);
    keys := q eq 2 select [1362434, 4089478963] else [188733, 34869438332];
    d := q eq 2 select 2 else 4;
    hs := [&+[R | v[j+1]*s^(d-j)*t^j : j in [0..d]]
        : v in VectorSpace(GF(q), d+1)];
    zs := q eq 2 select [&+[R | v[j+1]*s^(6-j)*t^j : j in [0..6]]
        : v in VectorSpace(GF(q), 7)] else [R!0];
    for key in keys do
        W := aInvariantsFromKey(R, key);
        canonical := orbit_invariant(W, invariant_context(W));
        least := q^(q eq 2 select 32 else 27);
        changed := false;
        for g in GL(2, GF(q)) do
            if q eq 3 and [c : c in Eltseq(g) | c ne 0][1] ne 1 then continue; end if;
            base := [g[1, 1]*s+g[1, 2]*t, g[2, 1]*s+g[2, 2]*t];
            A := [Evaluate(a, base) : a in W];
            for h in hs, z in zs do
                r := q eq 2 select h^2+A[1]*h else h;
                Z := equation_change(A, GF(q)!1, r, q eq 2 select h else R!0, z);
                forms := q eq 2 select [Z[1], Z[3], Z[4], Z[5]] else [Z[2], Z[4], Z[5]];
                degrees := q eq 2 select [2, 6, 8, 12] else [4, 8, 12];
                value := 0;
                for i in [1..#forms], j in [0..degrees[i]] do
                    value := q*value+Integers()!MonomialCoefficient(forms[i], s^(degrees[i]-j)*t^j);
                end for;
                least := Min(least, value);
                if not changed and FibrationKey(Z) ne key then
                    assert orbit_invariant(Z, invariant_context(Z)) eq canonical;
                    changed := true;
                end if;
            end for;
        end for;
        assert changed and least eq canonical;
    end for;
end procedure;
