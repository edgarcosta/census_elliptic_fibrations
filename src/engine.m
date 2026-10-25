function vector_numbers(M, p)
    weights := Matrix(Integers(),Ncols(M),1,
                      [ p^(Ncols(M)-i) : i in [1..Ncols(M)] ]);
    return Eltseq(ChangeRing(M,Integers())*weights);
end function;

function vector_number(v, p)
    assert IsPrime(p);
    n := 0;
    for a in Eltseq(v) do
        digit := Integers()!a;
        assert digit ge 0 and digit lt p;
        n := p*n+digit;
    end for;
    return n;
end function;

function orbit_representatives(G)
    // The vector of least key in each orbit of G on its natural module,
    // sorted by key.
    p := #BaseRing(G);
    reps := [];
    for O in Orbits(G) do
        m, i := Min([ vector_number(v,p) : v in O ]);
        Append(~reps, <m,O[i]>);
    end for;
    Sort(~reps, func< x,y | x[1]-y[1] >);
    return [ x[2] : x in reps ];
end function;

function affine_matrix(L, b)
    // [L,0; b,1], acting on rows (v,1) as v -> v*L+b.
    L := Matrix(L);
    k := BaseRing(L);
    d := Nrows(L);
    X := DiagonalJoin(L, IdentityMatrix(k,1));
    for j in [1..d] do X[d+1,j] := b[j]; end for;
    return GL(d+1,k)!X;
end function;

function affine_parts(X, d)
    X := Matrix(X);
    return Submatrix(X,1,1,d,d), Vector(BaseRing(X), [ X[d+1,j] : j in [1..d] ]);
end function;

function affine_quotient(E, S)
    // The shift subspace is stable under the linear parts. Projection and
    // lifting matrices act on rows; a lift is inverse only modulo shifts.
    n := Degree(E)-1;
    k := BaseRing(E);
    if Ngens(E) eq 0 then
        V := VectorSpace(k,n);
        Q, phi := quo< V | S >;
        d := Dimension(Q);
        return sub< GL(d+1,k) | >, d,
               Matrix(k, n, d, &cat[ Eltseq(phi(V.i)) : i in [1..n] ]),
               Matrix(k, d, n, &cat[ Eltseq(Q.i @@ phi) : i in [1..d] ]);
    end if;
    parts := [];
    for i in [1..Ngens(E)] do
        L, b := affine_parts(E.i,n);
        Append(~parts, <L,b>);
    end for;
    M := GModule(E, [ x[1] : x in parts ]);
    T := sub< M | [ M!Eltseq(v) : v in Basis(S) ] >;
    assert Dimension(T) eq Dimension(S);
    Q, proj := quo< M | T >;
    d := Dimension(Q);
    gens := [ affine_matrix(ActionGenerator(Q,i), Eltseq(proj(M!Eltseq(parts[i][2]))))
              : i in [1..Ngens(E)] ];
    toQ := Matrix(k, n, d, &cat[ Eltseq(proj(M.i)) : i in [1..n] ]);
    fromQ := Matrix(k, d, n, &cat[ Eltseq(Q.i @@ proj) : i in [1..d] ]);
    return sub< GL(d+1,k) | gens >, d, toQ, fromQ;
end function;

function project_subspaces(subspaces, toQ)
    // The distinct images of affine subspaces <x,U> under the matrix toQ.
    k := BaseRing(toQ);
    W := VectorSpace(k, Ncols(toQ));
    images := [* *];
    for xU in subspaces do
        x, U := Explode(xU);
        UQ := sub< W | [ W!(u*toQ) : u in Basis(U) ] >;
        pair := <ReduceVector(UQ,W!(x*toQ)),UQ>;
        if &and[ pair ne A : A in images ] then Append(~images,pair); end if;
    end for;
    return images;
end function;

FIBRE_FORMAT := recformat< Dimension, Group, Loci, Keys >;

function orbit_fibre(V, E, T, loci, keys)
    E,d,toQ,fromQ := affine_quotient(E,T);
    projected := [* project_subspaces(locus,toQ) : locus in loci *];
    return rec< FIBRE_FORMAT | Dimension := d, Group := E, Loci := projected,
        Keys := func< P | keys(P*fromQ) > >;
end function;

function affine_intersection(A, B)
    x,U := Explode(A);
    y,V := Explode(B);
    k := BaseRing(U);
    d := Degree(U);
    BU := Dimension(U) eq 0 select ZeroMatrix(k,0,d) else BasisMatrix(U);
    BV := Dimension(V) eq 0 select ZeroMatrix(k,0,d) else BasisMatrix(V);
    ok,c := IsConsistent(VerticalJoin(BU,-BV),y-x);
    if not ok then return false,_; end if;
    z := x+Vector(k,[ c[i] : i in [1..Dimension(U)] ])*BU;
    W := U meet V;
    return true,<ReduceVector(W,z),W>;
end function;

function union_intersection(A, B)
    result := [* *];
    for X in A do
        for Y in B do
            ok,Z := affine_intersection(X,Y);
            if ok and &and[ Z ne A : A in result ] then Append(~result,Z); end if;
        end for;
    end for;
    return result;
end function;

function affine_fixed_count(g, d, A)
    x,U := Explode(A);
    k := BaseRing(U);
    L,b := affine_parts(g,d);
    L -:= IdentityMatrix(k,d);
    rhs := -(x*L+b);
    if Dimension(U) eq 0 then return IsZero(rhs) select 1 else 0; end if;
    D := BasisMatrix(U)*L;
    return IsConsistent(D,rhs) select (#k)^(Dimension(U)-Rank(D)) else 0;
end function;

function subset_orbit_count(E, d, points)
    if #points eq 0 then return 0; end if;
    k := BaseRing(E);
    P := Matrix(k,#points,d,&cat[ Eltseq(v) : v in points ]);
    ones := Matrix(k,#points,1,[ k | 1 : i in [1..#points] ]);
    total := 0;
    for g in E do
        L,b := affine_parts(g,d);
        Y := P*(L-IdentityMatrix(k,d))+ones*Matrix(k,1,d,Eltseq(b));
        total +:= #[ i : i in [1..#points] | IsZero(Y[i]) ];
    end for;
    assert total mod #E eq 0;
    return total div #E;
end function;

procedure add_signed_term(~terms, A, n)
    for i in [1..#terms] do
        if terms[i][1] eq A then
            terms[i][2] +:= n;
            return;
        end if;
    end for;
    Append(~terms,<A,n>);
end procedure;

function union_orbit_count(E, d, components)
    if #components eq 0 then return 0; end if;
    if #components eq 1 then
        total := &+[ affine_fixed_count(g,d,components[1]) : g in E ];
    elif &+[ (#BaseRing(E))^Dimension(A[2]) : A in components ] le 1000000 then
        points := {};
        for A in components do points join:= { A[1]+v : v in A[2] }; end for;
        return subset_orbit_count(E,d,SetToSequence(points));
    else
        // Signed intersections correct overlaps without expanding large spaces.
        terms := [* *];
        for A in components do
            additions := [* <A,1> *];
            for term in terms do
                ok,C := affine_intersection(A,term[1]);
                if ok then add_signed_term(~additions,C,-term[2]); end if;
            end for;
            for term in additions do add_signed_term(~terms,term[1],term[2]); end for;
            terms := [* term : term in terms | term[2] ne 0 *];
        end for;
        total := &+[ Integers() | term[2]*&+[ affine_fixed_count(g,d,term[1]) : g in E ]
            : term in terms ];
    end if;
    assert total mod #E eq 0;
    return total div #E;
end function;

function membership_pattern(loci, v)
    return < &or[ (v-A[1]) in A[2] : A in locus ] : locus in loci >;
end function;

function fibre_counts(F)
    k := BaseRing(F`Group);
    W := VectorSpace(k,F`Dimension);
    m := #F`Loci;
    subsets := [ Setseq(S) : S in Subsets({1..m}) ];
    Sort(~subsets,func< A,B | #B-#A >);
    exact := AssociativeArray();
    for indices in subsets do
        locus := [* <W!0,W> *];
        for i in indices do locus := union_intersection(locus,F`Loci[i]); end for;
        n := union_orbit_count(F`Group,F`Dimension,locus);
        for pattern in Keys(exact) do
            if &and[ pattern[i] : i in indices ] then n -:= exact[pattern]; end if;
        end for;
        pattern := < i in indices : i in [1..m] >;
        assert n ge 0;
        exact[pattern] := n;
    end for;
    return exact;
end function;

function fibre_equation_counts(F)
    F`Group := sub< F`Group | >;
    return fibre_counts(F);
end function;

function fibre_representatives(F, selected)
    d := F`Dimension;
    E := F`Group;
    k := BaseRing(E);
    W := VectorSpace(k,d);
    result := [* *];
    batch := 4096;
    if (#k)^d gt 65536 then
        p := #k;
        parts := [* *];
        for g in E do
            if IsIdentity(g) then continue; end if;
            L,b := affine_parts(g,d);
            Append(~parts,<L,b>);
        end for;
        if #E gt 1 then seen := [ false : i in [1..p^d] ]; end if;
        powers := [ p^(d-i) : i in [1..d] ];
        for lo in [0..p^d-1 by batch] do
            hi := Min(lo+batch-1,p^d-1);
            P := Transpose(Matrix(k,d,hi-lo+1,
                &cat[ [ (n div power) mod p : n in [lo..hi] ] : power in powers ]));
            indices := [];
            if #E eq 1 then
                indices := [1..Nrows(P)];
            else
                ones := Matrix(k,hi-lo+1,1,[ k | 1 : n in [lo..hi] ]);
                images := [ vector_numbers(P*A[1]+ones*Matrix(A[2]),p) : A in parts ];
                for n in [lo..hi] do
                    if seen[n+1] then continue; end if;
                    Append(~indices,n-lo+1);
                    seen[n+1] := true;
                    for image in images do seen[image[n-lo+1]+1] := true; end for;
                end for;
            end if;
            patterns := [* membership_pattern(F`Loci,W!Eltseq(P[i])) : i in indices *];
            keep := [ i : i in [1..#indices] | patterns[i] in selected ];
            if #keep eq 0 then continue; end if;
            keys := F`Keys(RowSubmatrix(P,[ indices[i] : i in keep ]));
            for j in [1..#keep] do Append(~result,<keys[j],patterns[keep[j]]>); end for;
        end for;
        return result;
    end if;
    if #E eq 1 or d eq 0 then
        reps := [ v : v in W ];
    else
        parts := [* *];
        for i in [1..Ngens(E)] do
            L,b := affine_parts(E.i,d);
            Append(~parts,<L,b>);
        end for;
        ok,x0 := IsConsistent(HorizontalJoin([ A[1]-IdentityMatrix(k,d) : A in parts ]),
            Vector(k,&cat[ Eltseq(-A[2]) : A in parts ]));
        if ok then
            reps := [ x0+O[1] : O in Orbits(sub< GL(d,k) | [ A[1] : A in parts ] >) ];
        else
            reps := [ W![ O[1][j] : j in [1..d] ] : O in Orbits(E) | O[1][d+1] eq 1 ];
        end if;
    end if;
    for lo in [1..#reps by batch] do
        vectors := reps[lo..Min(lo+batch-1,#reps)];
        patterns := [* membership_pattern(F`Loci,v) : v in vectors *];
        keep := [ i : i in [1..#vectors] | patterns[i] in selected ];
        if #keep eq 0 then continue; end if;
        P := Matrix(k,#keep,d,&cat[ Eltseq(vectors[i]) : i in keep ]);
        keys := F`Keys(P);
        for j in [1..#keep] do Append(~result,<keys[j],patterns[keep[j]]>); end for;
    end for;
    return result;
end function;
