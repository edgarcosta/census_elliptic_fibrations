// Staged orbits of y^2 = x^3 + C*x^2 + A*x + B over GF(3).
// Load ef_common.m and ef3_models.m first. A model is the triple [C,A,B],
// with degrees
// 4,8,12.  Prefixes are tuples <C,A,H>, where H consists of pairs <M,q>
// that fix C,A under base substitution followed by x -> x+q.
//
// These routines are deliberately separate from the minimality filter.
// FibrationFiberOrbitKeysF3 is a reference/list interface.  A full census
// should stream one prefix at a time, without retaining all returned keys.

function EF3_QuarticOrbitData(R)
    k := BaseRing(R);
    assert #k eq 3 and Rank(R) eq 2;
    matrices := EF3_BaseMatrices(k);
    seen := [ false : i in [1..3^5] ];
    result := [* *];
    for n in [0..3^5-1] do
        if seen[n+1] then continue; end if;
        C := EF3_FormFromKey(R,4,n);
        H := [* *];
        for M in matrices do
            CM := EF_BaseChange(C,M);
            seen[EF3_Key(EF_Coefficients(CM,4))+1] := true;
            if CM eq C then Append(~H,M); end if;
        end for;
        Append(~result,<C,H>);
    end for;
    return result;
end function;

function EF3_DivideByQuartic(C,f)
    // Exact homogeneous division f/C, for deg(C)=4 and deg(f)=8.
    // A triangular convolution solve avoids any multivariate division
    // convention.  The final assertion checks all coefficients.
    R := Parent(C);
    k := BaseRing(R);
    assert C ne 0;
    cc := EF_Coefficients(C,4);
    ff := EF_Coefficients(f,8);
    first := Min([ i : i in [1..5] | cc[i] ne 0 ]);
    qq := [ k!0 : i in [1..5] ];
    for h in [0..4] do
        value := ff[first+h];
        for j in [0..h-1] do
            index := first+h-j;
            if index le 5 then
                value -:= cc[index]*qq[j+1];
            end if;
        end for;
        qq[h+1] := value/cc[first];
    end for;
    q := EF_FormFromVector(R,4,qq);
    assert C*q eq f;
    return q;
end function;

function FibrationOrbitPrefixesF3(R)
    // C!=0: quotient A by C*V4, a space of dimension 4, and then take
    // its stabilizer orbits.  C=0: take ordinary orbits of A in V8.
    // Expected counts: 933 prefixes with C!=0 and 927 with C=0.
    k := BaseRing(R);
    assert #k eq 3 and Rank(R) eq 2;
    V8 := VectorSpace(k,9);
    B4 := EF_BinaryBasis(R,4);
    prefixes := [* *];
    for item in EF3_QuarticOrbitData(R) do
        C := item[1];
        HC := item[2];
        if C eq 0 then
            seen := [ false : i in [1..3^9] ];
            for n in [0..3^9-1] do
                if seen[n+1] then continue; end if;
                A := EF3_FormFromKey(R,8,n);
                H := [* *];
                for M in HC do
                    AM := EF_BaseChange(A,M);
                    seen[EF3_Key(EF_Coefficients(AM,8))+1] := true;
                    if AM eq A then Append(~H,<M,R!0>); end if;
                end for;
                Append(~prefixes,<C,A,H>);
            end for;
        else
            S := sub< V8 |
                [ V8!EF_Coefficients(C*r,8) : r in B4 ] >;
            assert Dimension(S) eq 5;
            Q := sub< V8 |
                [ ReduceVector(S,V8.i) : i in [1..9] ] >;
            assert Dimension(Q) eq 4;
            QB := Basis(Q);
            W := VectorSpace(k,4);
            seen := [ false : i in [1..3^4] ];
            for n in [0..3^4-1] do
                if seen[n+1] then continue; end if;
                v := EF3_VectorFromKey(W,n);
                avec := V8!(&+[ v[i]*QB[i] : i in [1..4] ]);
                A := EF_FormFromVector(R,8,avec);
                H := [* *];
                for M in HC do
                    AM := EF_BaseChange(A,M);
                    w := ReduceVector(S,V8!EF_Coefficients(AM,8));
                    j := EF3_Key(Coordinates(Q,w));
                    seen[j+1] := true;
                    if w eq avec then
                        // AM+2*C*q=A, and 2^(-1)=2 in GF(3).
                        q := 2*EF3_DivideByQuartic(C,A-AM);
                        assert AM+2*C*q eq A;
                        Append(~H,<M,q>);
                    end if;
                end for;
                Append(~prefixes,<C,A,H>);
            end for;
        end if;
    end for;
    return prefixes;
end function;

function EF3_FiberData(prefix)
    // Return Q, maps, S, where Q complements the translation image S
    // in V12 and maps contains the distinct affine maps <L,b> on Q.
    // Coefficient vectors act on matrices from the left: v -> v*L+b.
    //
    // C!=0: no nonzero x-translation fixes the chosen A, so S=0.
    // C=0: the translations add r^3+A*r, which is GF(3)-linear in r.
    assert #prefix eq 3;
    C := prefix[1];
    A := prefix[2];
    H := prefix[3];
    R := Parent(C);
    k := BaseRing(R);
    assert #k eq 3;
    V12 := VectorSpace(k,13);
    if C eq 0 then
        S := sub< V12 |
            [ V12!EF_Coefficients(r^3+A*r,12)
              : r in EF_BinaryBasis(R,4) ] >;
        assert Dimension(S) in {4,5};
    else
        S := sub< V12 | [] >;
    end if;
    Q := sub< V12 |
        [ ReduceVector(S,V12.i) : i in [1..13] ] >;
    QB := Basis(Q);
    d := Dimension(Q);
    assert d eq 13-Dimension(S);
    assert d le 13;
    W := VectorSpace(k,d);
    maps := [* *];
    for h in H do
        M := h[1];
        q := h[2];
        assert EF_BaseChange(C,M) eq C;
        AM := EF_BaseChange(A,M);
        assert AM+2*C*q eq A;
        origin := q^3+AM*q+C*q^2;
        origin := ReduceVector(S,V12!EF_Coefficients(origin,12));
        shift := W!Coordinates(Q,origin);
        rows := [];
        for z in QB do
            B := EF_FormFromVector(R,12,z);
            BM := EF_BaseChange(B,M);
            w := ReduceVector(S,V12!EF_Coefficients(BM,12));
            Append(~rows,Coordinates(Q,w));
        end for;
        linear := Matrix(k,d,d,&cat rows);
        duplicate := false;
        for T in maps do
            if T[1] eq linear and T[2] eq shift then
                duplicate := true;
                break;
            end if;
        end for;
        if not duplicate then Append(~maps,<linear,shift>); end if;
    end for;
    assert #maps gt 0;
    return Q,maps,S;
end function;

function EF3_AffineBurnsideCount(maps)
    // Fixed points of v -> v*L+b solve v*(L-I)=-b.  Membership in the
    // row space decides solvability, and nullity gives the number.
    assert #maps gt 0;
    k := BaseRing(maps[1][1]);
    d := Nrows(maps[1][1]);
    W := VectorSpace(k,d);
    numerator := 0;
    for T in maps do
        L := T[1]-IdentityMatrix(k,d);
        image := sub< W |
            [ W![ L[i,j] : j in [1..d] ] : i in [1..d] ] >;
        if T[2] in image then
            numerator +:= (#k)^(d-Dimension(image));
        end if;
    end for;
    assert numerator mod #maps eq 0;
    return numerator div #maps;
end function;

function EF3_FiberBurnsideCount(prefix)
    Q,maps,S := EF3_FiberData(prefix);
    return EF3_AffineBurnsideCount(maps);
end function;

function FibrationStagedBurnsideCountF3(R : Prefixes := [])
    // All models, including Delta=0 and nonminimal models.  This checks
    // enumeration without traversing the 1.3 billion representative set.
    prefixes := Prefixes;
    if #prefixes eq 0 then prefixes := FibrationOrbitPrefixesF3(R); end if;
    count := 0;
    for prefix in prefixes do
        count +:= EF3_FiberBurnsideCount(prefix);
    end for;
    return count;
end function;

function SingularOrbitCountF3(R)
    // If C!=0 and Delta=0, the repeated root r=A/C is integral over the
    // normal ring GF(3)[s,t], hence is a quartic polynomial.  We have
    // A=C*r and B=C*r^2-r^3; x -> x+r sends the model to (C,0,0).
    // Thus this part has one orbit for each nonzero quartic C orbit.
    // If C=0, Delta=0 iff A=0, and B is taken modulo quartic cubes.
    quartics := EF3_QuarticOrbitData(R);
    nonzero := #quartics-1;
    H := [* *];
    for M in EF3_BaseMatrices(BaseRing(R)) do
        Append(~H,<M,R!0>);
    end for;
    zero := EF3_FiberBurnsideCount(<R!0,R!0,H>);
    return nonzero+zero,<nonzero,zero>;
end function;

function EF3_PrefixGenericSmooth(prefix,B)
    // Precondition: prefix was produced by FibrationOrbitPrefixesF3.
    // When C!=0, its A is reduced modulo C*V4.  The singular-root argument
    // above then gives Delta=0 iff A=B=0.  When C=0, Delta=-A^3.
    // This avoids constructing a discriminant for every representative.
    if prefix[1] eq 0 then return prefix[2] ne 0; end if;
    return prefix[2] ne 0 or B ne 0;
end function;

function FibrationPrefixStatisticsF3(prefixes)
    // Histogram keys are <C_is_zero, quotient_dimension, distinct_maps>.
    // The two point counts include, or skip, the entirely singular
    // prefix C=A=0. They count scanned vectors, not orbit representatives.
    histogram := AssociativeArray();
    allPoints := 0;
    smoothPoints := 0;
    for prefix in prefixes do
        Q,maps,S := EF3_FiberData(prefix);
        d := Dimension(Q);
        key := <prefix[1] eq 0,d,#maps>;
        found,count := IsDefined(histogram,key);
        if not found then count := 0; end if;
        histogram[key] := count+1;
        allPoints +:= 3^d;
        if prefix[1] ne 0 or prefix[2] ne 0 then
            smoothPoints +:= 3^d;
        end if;
    end for;
    return histogram,allPoints,smoothPoints;
end function;

function FibrationFiberOrbitKeysF3(prefix : GenericSmoothOnly := true)
    // One key per orbit, and <quotient dimension, distinct maps,
    // all orbits, retained orbits>.  A key has 27 base-3 digits C||A||B.
    // Precondition: prefix was produced by FibrationOrbitPrefixesF3.
    // This allocates a list for one prefix; use a streaming driver for
    // the complete census and do not concatenate all these key lists.
    C := prefix[1];
    A := prefix[2];
    R := Parent(C);
    k := BaseRing(R);
    if C ne 0 then
        V8 := VectorSpace(k,9);
        SC := sub< V8 | [ V8!EF_Coefficients(C*r,8)
                         : r in EF_BinaryBasis(R,4) ] >;
        avec := V8!EF_Coefficients(A,8);
        assert ReduceVector(SC,avec) eq avec;
    end if;
    Q,maps,S := EF3_FiberData(prefix);
    QB := Basis(Q);
    d := Dimension(Q);
    W := VectorSpace(k,d);
    V12 := VectorSpace(k,13);
    seen := [ false : i in [1..3^d] ];
    prefixDigits := EF_Coefficients(C,4) cat EF_Coefficients(A,8);
    keys := [ Integers() | ];
    total := 0;
    for n in [0..3^d-1] do
        if seen[n+1] then continue; end if;
        v := EF3_VectorFromKey(W,n);
        for T in maps do
            w := v*T[1]+T[2];
            seen[EF3_Key(Eltseq(w))+1] := true;
        end for;
        total +:= 1;
        z := V12!(&+[ v[i]*QB[i] : i in [1..d] ]);
        B := EF_FormFromVector(R,12,z);
        keep := true;
        if GenericSmoothOnly then
            keep := EF3_PrefixGenericSmooth(prefix,B);
        end if;
        if keep then
            Append(~keys,EF3_Key(prefixDigits cat Eltseq(z)));
        end if;
    end for;
    return keys,<d,#maps,total,#keys>;
end function;
