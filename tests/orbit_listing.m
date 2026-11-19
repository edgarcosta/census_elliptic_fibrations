import "tests/support/_weierstrass_change.m": homogeneous_forms, equation_change;

function fibre_coordinates(W)
    R := Universe(W);
    indices := Characteristic(R) eq 2 select [4, 5] else [5];
    degrees := [2, 4, 6, 8, 12];
    return Vector(BaseRing(R), &cat[[MonomialCoefficient(W[i],
        R.1^(degrees[i]-j)*R.2^j) : j in [0..degrees[i]]] : i in indices]);
end function;

function fibre_model(W, v)
    R := Universe(W);
    offset := 0;
    for i in (Characteristic(R) eq 2 select [4, 5] else [5]) do
        d := i eq 4 select 8 else 12;
        W[i] := &+[v[offset+j+1]*R.1^(d-j)*R.2^j : j in [0..d]];
        offset +:= d+1;
    end for;
    return W;
end function;

function independent_fibre_action(R, prefix)
    p := Characteristic(R);
    k := BaseRing(R);
    W := p eq 2 select [R | prefix[1], 0, prefix[2], 0, 0]
        else [R | 0, prefix[1], 0, prefix[2], 0];
    prefix_indices := p eq 2 select [1, 2, 3] else [1, 2, 3, 4];
    same_prefix := func<A | [A[j] : j in prefix_indices]
        eq [W[j] : j in prefix_indices]>;
    V := Parent(fibre_coordinates(W));
    images := {};
    for h in homogeneous_forms(R, p eq 2 select 6 else 4) do
        A := p eq 2 select equation_change(W, k!1, R!0, R!0, h)
            else equation_change(W, k!1, h, R!0, R!0);
        if same_prefix(A) then Include(~images, fibre_coordinates(A)); end if;
    end for;
    T := sub<V | Setseq(images)>;
    assert #images eq #T;
    Q, project := quo<V | T>;
    d := Dimension(Q);
    projection := Matrix(k, Dimension(V), d,
        &cat[Eltseq(project(v)) : v in Basis(V)]);
    lifts := [fibre_model(W, v @@ project) : v in Basis(Q)];
    matrices := {};
    // Units act trivially on the normalized coefficients over both fields.
    for g in GL(2, k) do
        base := [g[1, 1]*R.1+g[1, 2]*R.2, g[2, 1]*R.1+g[2, 2]*R.2];
        A := [Evaluate(a, base) : a in W];
        for h in homogeneous_forms(R, p eq 2 select 2 else 4) do
            r := p eq 2 select h^2+A[1]*h else h;
            q := p eq 2 select h else R!0;
            changed := equation_change(A, k!1, r, q, R!0);
            if not same_prefix(changed) then continue; end if;
            b := project(fibre_coordinates(changed));
            rows := [];
            for lift in lifts do
                B := [Evaluate(a, base) : a in lift];
                row := project(fibre_coordinates(equation_change(B, k!1, r, q, R!0)))-b;
                Append(~rows, Eltseq(row) cat [k!0]);
            end for;
            Append(~rows, Eltseq(b) cat [k!1]);
            Include(~matrices, GL(d+1, k)!Matrix(rows));
            // Other restoring quartics differ by a translation already in T.
            if p eq 3 then break; end if;
        end for;
    end for;
    G := sub<GL(d+1, k) | Setseq(matrices)>;
    assert #G eq #matrices;
    return G, projection;
end function;

// Disjoint representative orbits exhaust the independently constructed fibre.
procedure test_orbit_listing(p, indices)
    R<s, t> := PolynomialRing(GF(p), 2);
    prefixes := FibrationPrefixes(R);
    for i in indices do
        G, projection := independent_fibre_action(R, prefixes[i]);
        n := Nrows(projection);
        d := Ncols(projection);
        keys := FibrationRepresentativeKeys(R, i);
        assert #keys gt 0;
        degrees := p eq 2 select [2, 6] else [4, 8];
        prefix_key := 0;
        for j in [1..2], a in [0..degrees[j]] do
            prefix_key := p*prefix_key+Integers()!MonomialCoefficient(
                prefixes[i][j], s^(degrees[j]-a)*t^a);
        end for;
        offset := prefix_key*p^n;
        assert &and[offset le key and key lt offset+p^n : key in keys];
        weights := Matrix(Integers(), d+1, 1, [p^(d-j) : j in [1..d]] cat [0]);
        seen := [false : j in [1..p^d]];
        covered := 0;
        for lo in [1..#keys by 4096] do
            batch := keys[lo..Min(lo+4095, #keys)];
            P := Matrix(GF(p), #batch, n,
                &cat[[(key div p^(n-j)) mod p : j in [1..n]] : key in batch]);
            P := HorizontalJoin(P*projection, Matrix(GF(p), #batch, 1, [1 : key in batch]));
            orbit_keys := [Eltseq(ChangeRing(P*Matrix(g), Integers())*weights) : g in G];
            for j in [1..#batch] do
                orbit := {image[j]+1 : image in orbit_keys};
                error if &or[seen[a] : a in orbit],
                    Sprintf("GF(%o), prefix %o: repeated representative orbit at key %o", p, i, batch[j]);
                for a in orbit do seen[a] := true; end for;
                covered +:= #orbit;
            end for;
        end for;
        error if covered ne p^d,
            Sprintf("GF(%o), prefix %o: representative orbits cover %o of %o fibre points", p, i, covered, p^d);
    end for;
end procedure;

test_orbit_listing(2, [1, 2, 88]);
test_orbit_listing(3, [1, 288, 1330]);
