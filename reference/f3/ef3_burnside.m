// Independent Burnside count for characteristic-3 Weierstrass models.
// No other project files are required.
//
// Models are triples [C,A,B] in GF(3)[s,t], of degrees 4,8,12:
//     y^2 = x^3 + C*x^2 + A*x + B.
// A group triple <M,u,q> applies base pullback f(z) -> f(M*z), then
// coefficient scaling (C,A,B) -> (u^2*C,u^4*A,u^6*B), then x -> x+q.
// Thus, with C0=u^2*C(M*z), A0=u^4*A(M*z), B0=u^6*B(M*z),
//     C' = C0;
//     A' = A0 + 2*C0*q;
//     B' = B0 + A0*q + C0*q^2 + q^3.
//
// Right action convention: apply g, then h, to obtain g*h. Consequently
//   <M,u,q>*<N,v,r> = <M*N,u*v,v^2*q(N*z)+r>.
// The inverse is <M^-1,u^-1,-u^-2*q(M^-1*z)>.
// Over GF(3), u=+/-1 acts trivially on all coefficients. Also M and -M
// have identical actions. The effective group is V4 semidirect PGL(2,3),
// of order 3^5*24=5832. Its action on the 27 coefficient coordinates is
// affine. Appending a final coordinate 1 realizes it in GL(28,3), acting
// on row vectors by right multiplication.
//
// Each coefficient block is ordered s^d,s^(d-1)*t,...,t^d.
// If w' = w*L+b, the 28 by 28 matrix is [ L  0 ; b  1 ].

function EFB_BinaryBasis(R, d)
    return [ R | R.1^(d-i)*R.2^i : i in [0..d] ];
end function;

function EFB_CheckRing(R)
    assert Type(R) eq RngMPol;
    assert Rank(R) eq 2;
    k := BaseRing(R);
    assert Type(k) eq FldFin;
    assert #k eq 3;
    return k;
end function;

function EFB_Coefficients(f, d)
    R := Parent(f);
    return [ BaseRing(R) | MonomialCoefficient(f,m)
             : m in EFB_BinaryBasis(R,d) ];
end function;

function EFB_Pullback(f, M)
    R := Parent(f);
    s := R.1; t := R.2;
    return Evaluate(f,
        [ M[1,1]*s+M[1,2]*t, M[2,1]*s+M[2,2]*t ]);
end function;

function EFB_CheckGroupTriple(g)
    assert #g eq 3;
    M,u,q := Explode(g);
    R := Parent(q);
    k := EFB_CheckRing(R);
    assert Nrows(M) eq 2 and Ncols(M) eq 2;
    assert BaseRing(M) eq k;
    assert Determinant(M) ne 0;
    u := k!u;
    assert u ne 0;
    for m in Monomials(q) do
        assert TotalDegree(m) eq 4;
    end for;
    return R,k;
end function;

function EFB_Compose(g, h)
    // H(g*h)=H(g)*H(h): g is applied first.
    R,k := EFB_CheckGroupTriple(g);
    R1,k1 := EFB_CheckGroupTriple(h);
    assert R1 eq R;
    M,u,q := Explode(g);
    N,v,r := Explode(h);
    return <M*N, (k!u)*(k!v), (k!v)^2*EFB_Pullback(q,N)+r>;
end function;

function EFB_Inverse(g)
    R,k := EFB_CheckGroupTriple(g);
    M,u,q := Explode(g);
    N := M^-1;
    return <N, (k!u)^-1, -(k!u)^-2*EFB_Pullback(q,N)>;
end function;

function EFB_CanonicalBaseMatrix(M)
    // Over GF(3), choose the representative of {M,-M} whose first
    // nonzero entry in row order is 1.
    assert #BaseRing(M) eq 3;
    assert Determinant(M) ne 0;
    for c in Eltseq(M) do
        if c ne 0 then return c^-1*M; end if;
    end for;
    error "The zero matrix is not invertible";
end function;

function EFB_CanonicalTriple(g)
    R,k := EFB_CheckGroupTriple(g);
    M,u,q := Explode(g);
    // Removing u and changing M to -M do not change q over GF(3).
    return <EFB_CanonicalBaseMatrix(M), k!1, q>;
end function;

function EFB_Transform(CAB, g)
    R,k := EFB_CheckGroupTriple(g);
    assert #CAB eq 3;
    degrees := [4,8,12];
    for i in [1..3] do
        assert Parent(CAB[i]) eq R;
        for m in Monomials(CAB[i]) do
            assert TotalDegree(m) eq degrees[i];
        end for;
    end for;
    M,u,q := Explode(g);
    u := k!u;
    C := u^2*EFB_Pullback(CAB[1],M);
    A := u^4*EFB_Pullback(CAB[2],M);
    B := u^6*EFB_Pullback(CAB[3],M);
    return [ R | C, A+2*C*q, B+A*q+C*q^2+q^3 ];
end function;

function EFB_ModelRow(CAB)
    assert #CAB eq 3;
    k := EFB_CheckRing(Parent(CAB[1]));
    return Vector(k, EFB_Coefficients(CAB[1],4)
                     cat EFB_Coefficients(CAB[2],8)
                     cat EFB_Coefficients(CAB[3],12) cat [k!1]);
end function;

function EFB_AffineMatrix(g)
    // Build rows from transformed monomials, independently of any
    // prefix normalization or quotient-space code.
    R,k := EFB_CheckGroupTriple(g);
    M,u,q := Explode(g);
    u := k!u;
    rows := [];
    for m in EFB_BinaryBasis(R,4) do
        f := u^2*EFB_Pullback(m,M);
        row := EFB_Coefficients(f,4)
               cat EFB_Coefficients(2*f*q,8)
               cat EFB_Coefficients(f*q^2,12) cat [k!0];
        Append(~rows,row);
    end for;
    for m in EFB_BinaryBasis(R,8) do
        f := u^4*EFB_Pullback(m,M);
        row := [k!0 : i in [1..5]] cat EFB_Coefficients(f,8)
               cat EFB_Coefficients(f*q,12) cat [k!0];
        Append(~rows,row);
    end for;
    for m in EFB_BinaryBasis(R,12) do
        f := u^6*EFB_Pullback(m,M);
        row := [k!0 : i in [1..14]] cat EFB_Coefficients(f,12)
               cat [k!0];
        Append(~rows,row);
    end for;
    Append(~rows, [k!0 : i in [1..14]]
                  cat EFB_Coefficients(q^3,12) cat [k!1]);
    return Matrix(k,28,28,&cat rows);
end function;

function EFB_BaseRepresentatives(R)
    k := EFB_CheckRing(R);
    reps := [];
    for g in GL(2,k) do
        M := Matrix(g);
        if M eq EFB_CanonicalBaseMatrix(M) then
            Append(~reps,M);
        end if;
    end for;
    assert #reps eq 24;
    return reps;
end function;

function EFB_Quartics(R)
    k := EFB_CheckRing(R);
    B4 := EFB_BinaryBasis(R,4);
    return [ R | &+[ R | v[i]*B4[i] : i in [1..5] ]
             : v in VectorSpace(k,5) ];
end function;

function EFB_AffineGroupF3(R)
    k := EFB_CheckRing(R);
    I2 := IdentityMatrix(k,2);
    S := Matrix(k,2,2,[0,-1,1,0]);
    T := Matrix(k,2,2,[1,1,0,1]);
    D := Matrix(k,2,2,[-1,0,0,1]);
    triples := [* <S,k!1,R!0>, <T,k!1,R!0>, <D,k!1,R!0> *];
    for q in EFB_BinaryBasis(R,4) do
        Append(~triples,<I2,k!1,q>);
    end for;
    ambient := GL(28,k);
    images := [ ambient!EFB_AffineMatrix(g) : g in triples ];
    G := sub< ambient | images >;
    // These are actual GrpMatElt objects, not informal affine pairs.
    return G, [ G!h : h in images ], triples;
end function;

function EFB_FixedDimension(H)
    // Return -1 for an empty affine fixed locus, otherwise its dimension.
    // For row vectors, w*(L-I)=-b is solvable exactly when b belongs to
    // the row space of L-I. The sign does not affect this membership.
    k := BaseRing(H);
    assert Nrows(H) eq 28 and Ncols(H) eq 28;
    for i in [1..27] do assert H[i,28] eq 0; end for;
    assert H[28,28] eq 1;
    A := Submatrix(H,1,1,27,27)-IdentityMatrix(k,27);
    b := Matrix(k,1,27,[ H[28,j] : j in [1..27] ]);
    r := Rank(A);
    if Rank(VerticalJoin(A,b)) ne r then return -1; end if;
    return 27-r;
end function;

procedure EFB_Accumulate(H, ~spectrum, ~empty, ~fixedSum, ~elements)
    d := EFB_FixedDimension(H);
    if d eq -1 then
        empty +:= 1;
    else
        spectrum[d] +:= 1;
        fixedSum +:= 3^d;
    end if;
    elements +:= 1;
end procedure;

function EFB_BurnsideF3(R : UseGroupIterator := false,
                          CheckExpected := true)
    // Default: explicitly enumerate 24 projective matrices and 243
    // quartics, as in the original independent computation. This is
    // 5832 fixed-locus calculations on 27 variables; no model enumeration.
    // UseGroupIterator=true instead enumerates the generated matrix group.
    k := EFB_CheckRing(R);
    G,images,triples := EFB_AffineGroupF3(R);
    assert Order(G) eq 5832;
    spectrum := AssociativeArray(Integers());
    for d in [0..27] do spectrum[d] := 0; end for;
    empty := 0; fixedSum := 0; elements := 0;
    if UseGroupIterator then
        for h in G do
            EFB_Accumulate(Matrix(h),~spectrum,~empty,~fixedSum,~elements);
        end for;
    else
        I2 := IdentityMatrix(k,2);
        bases := [ EFB_AffineMatrix(<M,k!1,R!0>)
                   : M in EFB_BaseRepresentatives(R) ];
        shifts := [ EFB_AffineMatrix(<I2,k!1,q>) : q in EFB_Quartics(R) ];
        for H0 in bases do
            for Hq in shifts do
                EFB_Accumulate(H0*Hq,~spectrum,~empty,~fixedSum,~elements);
            end for;
        end for;
    end if;
    assert elements eq 5832;
    assert empty+(&+[ spectrum[d] : d in [0..27] ]) eq elements;
    assert fixedSum mod elements eq 0;
    orbits := fixedSum div elements;
    if CheckExpected then
        expected := AssociativeArray(Integers());
        for d in [0..27] do expected[d] := 0; end for;
        for pair in [ <27,1>, <15,81>, <13,242>, <10,216>,
                      <9,432>, <7,2592>, <5,432>, <3,972> ] do
            expected[pair[1]] := pair[2];
        end for;
        for d in [0..27] do assert spectrum[d] eq expected[d]; end for;
        assert empty eq 864;
        assert fixedSum eq 7627172630184;
        assert orbits eq 1307814237;
    end if;
    return orbits,spectrum,empty,fixedSum,G;
end function;

function EFB_PrefixBurnsideF3(R : CheckExpected := true)
    // Independent Burnside checks on the linear prefix representations.
    // Return <number of quartic C orbits, number of prefixes with C!=0,
    //         number of prefixes with C=0>.
    // The expected answer is <24,933,927>; the last two entries sum to
    // 1860, the number of orbits of all pairs (C,A) under the full group.
    // This function uses no staged prefix lists or normal-form choices.
    k := EFB_CheckRing(R);
    I2 := IdentityMatrix(k,2);
    I5 := IdentityMatrix(k,5);
    I9 := IdentityMatrix(k,9);
    I14 := IdentityMatrix(k,14);
    bases := [ EFB_AffineMatrix(<M,k!1,R!0>)
               : M in EFB_BaseRepresentatives(R) ];
    shifts := [ EFB_AffineMatrix(<I2,k!1,q>) : q in EFB_Quartics(R) ];
    cFixedSum := 0; aFixedSum := 0; prefixFixedSum := 0;
    for H0 in bases do
        L5 := Submatrix(H0,1,1,5,5);
        L9 := Submatrix(H0,6,6,9,9);
        cFixedSum +:= 3^(5-Rank(L5-I5));
        // On C=0, translations fix A, so only the 24 base maps act.
        aFixedSum +:= 3^(9-Rank(L9-I9));
        for Hq in shifts do
            L14 := Submatrix(H0*Hq,1,1,14,14);
            prefixFixedSum +:= 3^(14-Rank(L14-I14));
        end for;
    end for;
    assert cFixedSum mod 24 eq 0;
    assert aFixedSum mod 24 eq 0;
    assert prefixFixedSum mod 5832 eq 0;
    cOrbits := cFixedSum div 24;
    zeroPrefixes := aFixedSum div 24;
    allPrefixes := prefixFixedSum div 5832;
    nonzeroPrefixes := allPrefixes-zeroPrefixes;
    if CheckExpected then
        assert cOrbits eq 24;
        assert zeroPrefixes eq 927;
        assert allPrefixes eq 1860;
        assert nonzeroPrefixes eq 933;
    end if;
    return <cOrbits,nonzeroPrefixes,zeroPrefixes>;
end function;

procedure EFB_CheckActionConventions(R)
    // These checks specifically detect transposition, composition-order,
    // translation, and scalar-kernel errors.
    k := EFB_CheckRing(R);
    G,images,triples := EFB_AffineGroupF3(R);
    assert Order(G) eq 5832;
    baseGroup := sub< G | [images[i] : i in [1..3]] >;
    shiftGroup := sub< G | [images[i] : i in [4..8]] >;
    assert Order(baseGroup) eq 24;
    assert Order(shiftGroup) eq 243;
    assert IsAbelian(shiftGroup);
    assert IsNormal(G,shiftGroup);
    I28 := IdentityMatrix(k,28);
    for i in [1..#triples] do
        g := triples[i];
        assert EFB_AffineMatrix(EFB_Inverse(g)) eq Matrix(images[i])^-1;
        assert EFB_AffineMatrix(EFB_Compose(g,EFB_Inverse(g))) eq I28;
        for j in [1..#triples] do
            assert EFB_AffineMatrix(EFB_Compose(g,triples[j]))
                   eq Matrix(images[i])*Matrix(images[j]);
        end for;
    end for;
    s := R.1; t := R.2;
    sample := [ R | s^4+s^3*t+2*t^4,
                    s^8+2*s^5*t^3+s*t^7,
                    2*s^12+s^7*t^5+t^12 ];
    for g in triples do
        assert EFB_ModelRow(EFB_Transform(sample,g))
               eq EFB_ModelRow(sample)*EFB_AffineMatrix(g);
    end for;
    for lambda in k do
        if lambda eq 0 then continue; end if;
        for u in k do
            if u eq 0 then continue; end if;
            // In GF(3), these are exactly all four kernel triples.
            assert EFB_AffineMatrix(<lambda*IdentityMatrix(k,2),u,R!0>)
                   eq I28;
        end for;
    end for;
end procedure;
