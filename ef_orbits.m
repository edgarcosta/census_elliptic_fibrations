// Orbits of finite affine groups on fibres, shared by the census code for
// each characteristic. Load ef_common.m first. A fibre is an EF_Fibre:
//   Dimension:  d;
//   Group:      E <= GL(d+1,p) of matrices [L,0; b,1], acting on rows (v,1)
//               as v -> v*L+b;
//   Singular:   false, or <x0,U>, the E-invariant affine subspace x0+U of
//               fibre vectors whose model has Delta = 0;
//   NonMinimal: a sequence of affine subspaces <x,U> of the fibre whose
//               union is E-invariant and is the set of models that are not
//               minimal at some place (it may meet Singular);
//   Key:        a function from fibre vectors to model keys;
//   FromKey:    its inverse, from the key of any model over the prefix to
//               its fibre vector.

EF_Fibre := recformat< Dimension, Group, Singular, NonMinimal, Key, FromKey >;

function EF_FormAction(R, G, d)
    assert BaseRing(G) eq BaseRing(R);
    V := SymmetricPower(GModule(G),d);
    return hom< G -> GL(d+1,BaseRing(R)) |
        ActionGenerators(V) >;
end function;

function EF_OrbitRepresentatives(G)
    // The vector of least key in each orbit of G on its natural module,
    // sorted by key.
    p := #BaseRing(G);
    reps := [];
    for O in Orbits(G) do
        m, i := Min([ EF_PrimeFieldKey(v,p) : v in O ]);
        Append(~reps, <m,O[i]>);
    end for;
    Sort(~reps, func< x,y | x[1]-y[1] >);
    return [ x[2] : x in reps ];
end function;

function EF_AffineMatrix(L, b)
    // [L,0; b,1], acting on rows (v,1) as v -> v*L+b.
    L := Matrix(L);
    k := BaseRing(L);
    d := Nrows(L);
    X := DiagonalJoin(L, IdentityMatrix(k,1));
    for j in [1..d] do X[d+1,j] := b[j]; end for;
    return GL(d+1,k)!X;
end function;

function EF_AffineParts(X, d)
    X := Matrix(X);
    return Submatrix(X,1,1,d,d), Vector(BaseRing(X), [ X[d+1,j] : j in [1..d] ]);
end function;

function EF_AffineQuotient(E, S)
    // E <= GL(n+1,p) affine, S a subspace of F_p^n preserved by the linear
    // parts of E. Returns the affine group induced on the quotient module
    // F_p^n/S, its dimension d, the n x d matrix of the projection and the
    // d x n matrix of a lift, both acting on row vectors.
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
        L, b := EF_AffineParts(E.i,n);
        Append(~parts, <L,b>);
    end for;
    M := GModule(E, [ x[1] : x in parts ]);
    T := sub< M | [ M!Eltseq(v) : v in Basis(S) ] >;
    assert Dimension(T) eq Dimension(S);
    Q, proj := quo< M | T >;
    d := Dimension(Q);
    gens := [ EF_AffineMatrix(ActionGenerator(Q,i), Eltseq(proj(M!Eltseq(parts[i][2]))))
              : i in [1..Ngens(E)] ];
    toQ := Matrix(k, n, d, &cat[ Eltseq(proj(M.i)) : i in [1..n] ]);
    fromQ := Matrix(k, d, n, &cat[ Eltseq(Q.i @@ proj) : i in [1..d] ]);
    return sub< GL(d+1,k) | gens >, d, toQ, fromQ;
end function;

function EF_ProjectSubspaces(subspaces, toQ)
    // The distinct images of affine subspaces <x,U> under the matrix toQ.
    k := BaseRing(toQ);
    W := VectorSpace(k, Ncols(toQ));
    images := {};
    for xU in subspaces do
        x, U := Explode(xU);
        UQ := sub< W | [ W!(u*toQ) : u in Basis(U) ] >;
        Include(~images, <ReduceVector(UQ, W!(x*toQ)), UQ>);
    end for;
    return SetToSequence(images);
end function;

function EF_AffineFixedCount(X, d)
    L, b := EF_AffineParts(X,d);
    L -:= IdentityMatrix(BaseRing(L),d);
    if not IsConsistent(L,-b) then return 0; end if;
    return (#BaseRing(L))^(d-Rank(L));
end function;

function EF_AffineOrbitCount(E, d)
    total := &+[ EF_AffineFixedCount(g,d) : g in E ];
    assert total mod #E eq 0;
    return total div #E;
end function;

function EF_InSingular(F, v)
    if Type(F`Singular) eq BoolElt then return false; end if;
    x0, U := Explode(F`Singular);
    return (v-x0) in U;
end function;

function EF_SingularOrbitCount(F)
    if Type(F`Singular) eq BoolElt then return 0; end if;
    d := F`Dimension;
    x0, U := Explode(F`Singular);
    k := BaseRing(U);
    BU := Dimension(U) eq 0 select ZeroMatrix(k,0,d) else BasisMatrix(U);
    total := 0;
    for g in F`Group do
        L, b := EF_AffineParts(g,d);
        L -:= IdentityMatrix(k,d);
        // Fixed points x0+c*BU solve c*(BU*L) = -(x0*L+b).
        rhs := -(x0*L+b);
        if Dimension(U) eq 0 then
            if IsZero(rhs) then total +:= 1; end if;
        elif IsConsistent(BU*L, rhs) then
            total +:= (#k)^(Dimension(U)-Rank(BU*L));
        end if;
    end for;
    assert total mod #F`Group eq 0;
    return total div #F`Group;
end function;

function EF_SubsetOrbitCount(E, d, points)
    // Orbits of E on a finite E-invariant set of fibre vectors.
    if #points eq 0 then return 0; end if;
    k := BaseRing(E);
    n := #points;
    P := Matrix(k, n, d, &cat[ Eltseq(v) : v in points ]);
    ones := Matrix(k, n, 1, [ k | 1 : i in [1..n] ]);
    total := 0;
    for g in E do
        L, b := EF_AffineParts(g,d);
        Y := P*(L-IdentityMatrix(k,d)) + ones*Matrix(k,1,d,Eltseq(b));
        total +:= #[ i : i in [1..n] | IsZero(Y[i]) ];
    end for;
    assert total mod #E eq 0;
    return total div #E;
end function;

function EF_NonMinimalPoints(F)
    points := {};
    for xU in F`NonMinimal do
        x, U := Explode(xU);
        points join:= { x+u : u in U };
    end for;
    return SetToSequence(points);
end function;

function EF_NonMinimalNonsingular(F)
    points := EF_NonMinimalPoints(F);
    if Type(F`Singular) eq BoolElt or #points eq 0 then return points; end if;
    x0, U := Explode(F`Singular);
    k := BaseRing(U);
    d := F`Dimension;
    // v is in x0+U exactly when (v-x0)*N = 0, N spanning the annihilator of U.
    N := Transpose(BasisMatrix(Kernel(Transpose(
        Dimension(U) eq 0 select ZeroMatrix(k,1,d) else BasisMatrix(U)))));
    if Ncols(N) eq 0 then return []; end if;
    P := Matrix(k, #points, d, &cat[ Eltseq(v) : v in points ]);
    Y := P*N - Matrix(k, #points, 1, [ k | 1 : i in [1..#points] ])*Matrix(x0*N);
    return [ points[i] : i in [1..#points] | not IsZero(Y[i]) ];
end function;

function EF_FibreCounts(F)
    // <all orbits, orbits with Delta != 0, elliptic K3 orbits>.
    d := F`Dimension;
    E := F`Group;
    all := EF_AffineOrbitCount(E,d);
    smooth := all-EF_SingularOrbitCount(F);
    return <all, smooth, smooth-EF_SubsetOrbitCount(E,d,EF_NonMinimalNonsingular(F))>;
end function;

function EF_FibreRepresentatives(F)
    // Model keys of one representative per elliptic K3 orbit.
    d := F`Dimension;
    E := F`Group;
    k := BaseRing(E);
    W := VectorSpace(k,d);
    if Type(F`Singular) ne BoolElt and Dimension(F`Singular[2]) eq d then
        return [ Integers() | ];
    end if;
    bad := Set(EF_NonMinimalNonsingular(F));
    if #E eq 1 then
        reps := [ v : v in W ];
    else
        parts := [];
        for i in [1..Ngens(E)] do
            L, b := EF_AffineParts(E.i,d);
            Append(~parts, <L,b>);
        end for;
        ok, x0 := IsConsistent(
            HorizontalJoin([ x[1]-IdentityMatrix(k,d) : x in parts ]),
            Vector(k, &cat[ Eltseq(-x[2]) : x in parts ]));
        if ok then
            // Translating a fixed point to 0 makes the action linear.
            reps := [ x0+O[1] : O in Orbits(sub< GL(d,k) | [ x[1] : x in parts ] >) ];
        else
            // E preserves the last coordinate; keep the orbits where it is 1.
            reps := [ W![ O[1][j] : j in [1..d] ] : O in Orbits(E)
                      | O[1][d+1] eq 1 ];
        end if;
    end if;
    return [ F`Key(v) : v in reps | v notin bad and not EF_InSingular(F,v) ];
end function;

procedure EF_WriteRepresentatives(F, filename)
    // Appends 8-byte little-endian model keys, as ef3_stream.m writes them.
    keys := EF_FibreRepresentatives(F);
    file := Open(filename, "ab");
    for i in [1..#keys by 4096] do
        WriteBytes(file, &cat[ [ (key div 256^j) mod 256 : j in [0..7] ]
                               : key in keys[i..Min(i+4095,#keys)] ]);
    end for;
    Flush(file);
    delete file;
end procedure;

function EF_CensusCounts(prefixes, FibreOf)
    // Per-prefix counts and the totals <all, Delta != 0, elliptic K3>.
    counts := [ EF_FibreCounts(FibreOf(P)) : P in prefixes ];
    return counts, < &+[ c[1] : c in counts ], &+[ c[2] : c in counts ],
                     &+[ c[3] : c in counts ] >;
end function;
