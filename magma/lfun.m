LFUN_FORMAT := recformat<Polynomial, Degree, Sign, Lower, TraceSums,
    Height, TwistHeight, TamagawaExponent>;

function fibre_trace(a, delta, table, ordinary)
    a1, a2, a3 := Explode(a[1..3]);
    K := table`Field;
    if Characteristic(K) eq 2 then
        if a1 ne 0 then
            sign := IsZero(Trace(a2/a1^2+a3/a1^3)) select 1 else -1;
            return sign*ordinary[delta/a1^12];
        end if;
    else
        c := a2+a1^2;
        if c ne 0 then
            sign := IsSquare(c) select 1 else -1;
            return sign*ordinary[-delta/c^6];
        end if;
    end if;
    // Supersingular fibres are sparse on an ordinary generic curve.
    // Magma counts them directly, including everywhere supersingular inputs.
    return #K+1-#EllipticCurve(a);
end function;

procedure add_euler_sums(~sums, trace, norm, degree, good)
    previous := 2;
    current := trace;
    for r in [1..#sums div degree] do
        sums[r*degree] +:= degree*current;
        next := trace*current-(good select norm*previous else 0);
        previous := current;
        current := next;
    end for;
end procedure;

function euler_sums(a, delta, rows, p, limit)
    sums := [Integers() | 0 : n in [1..limit]];
    for row in rows do
        add_euler_sums(~sums, row`Trace, p^row`Degree, row`Degree, false);
    end for;
    for d in [1..limit] do
        table := FibrationTraceTable(p, d);
        K := table`Field;
        alpha := table`Generator;
        ordinary := AssociativeArray(K);
        power := K!1;
        representatives := [K |];
        order := #K-1;
        for e in [0..order-1] do
            ordinary[power] := table`Ordinary[e+1];
            orbit := [Integers() |];
            f := e;
            repeat
                Append(~orbit, f);
                f := (p*f) mod order;
            until f eq e;
            if #orbit eq d and e eq Min(orbit) then
                Append(~representatives, power);
            end if;
            power *:= alpha;
        end for;
        if d eq 1 then Append(~representatives, K!0); end if;
        for point in representatives do
            discriminant := Evaluate(delta, point);
            if discriminant eq 0 then continue; end if;
            values := [Evaluate(b, point) : b in a];
            trace := fibre_trace(values, discriminant, table, ordinary);
            add_euler_sums(~sums, trace, p^d, d, true);
        end for;
        if d eq 1 and not &or[row`Place cmpeq "inf" : row in rows] then
            weights := [2, 4, 6, 8, 12];
            values := [K!Coefficient(a[i], weights[i]) : i in [1..5]];
            trace := fibre_trace(values, K!Coefficient(delta, 24), table, ordinary);
            add_euler_sums(~sums, trace, p, 1, true);
        end if;
    end for;
    return sums;
end function;

intrinsic FibrationLPolynomial(E::CrvEll[FldFunRat]) -> Rec
{The signed L-polynomial record of an everywhere minimal elliptic K3 polynomial
model over GF(2)(t) or GF(3)(t), with Degree, Sign and Lower coefficients
c1 through c_floor(Degree/2). See also FibrationLocalData, FibrationTraceTable}
    F := BaseRing(E);
    k := BaseRing(F);
    require Type(k) eq FldFin and #k in [2, 3]:
        "constant field must be GF(2) or GF(3)";
    a := aInvariants(E);
    weights := [2, 4, 6, 8, 12];
    require &and[Denominator(a[i]) eq 1 and Degree(Numerator(a[i])) le weights[i]
        : i in [1..5]]:
        "coefficients must be polynomials of degrees at most 2,4,6,8,12";
    R<s, u> := PolynomialRing(k, 2);
    W := [R | &+[R | Coefficient(Numerator(a[i]), j)*s^j*u^(weights[i]-j)
        : j in [0..Degree(Numerator(a[i]))]] : i in [1..5]];
    kind, minimality := FibrationProperties(W);
    require kind eq "smooth" and minimality eq "everywhere minimal":
        "equation must define an everywhere minimal elliptic K3";
    p := #k;
    rows := FibrationLocalData(E);
    twists := FibrationLocalData(E : FrobeniusTwist := true);
    discriminant_degree := &+[r`Degree*r`DiscriminantValuation : r in rows];
    twist_degree := &+[r`Degree*r`DiscriminantValuation : r in twists];
    assert discriminant_degree eq 24 and twist_degree mod 12 eq 0;
    height := discriminant_degree div 12;
    twist_height := twist_degree div 12;
    tamagawa_exponent := &+[Valuation(r`Tamagawa, p) : r in twists]
        - &+[Valuation(r`Tamagawa, p) : r in rows];
    sign := IsEven(twist_height-height+tamagawa_exponent) select 1 else -1;
    degree := &+[r`Degree*r`Conductor : r in rows]-4;
    assert degree ge 0 and degree le 20;
    P := PolynomialRing(k);
    sums := euler_sums([P!b : b in a], P!Discriminant(E), rows, p, degree div 2);
    lower := [Integers() | 1];
    for n in [1..#sums] do
        numerator := &+[sums[j]*lower[n-j+1] : j in [1..n]];
        assert numerator mod n eq 0;
        Append(~lower, numerator div n);
    end for;
    coefficients := [Integers() | 0 : i in [0..degree]];
    for i in [0..degree div 2] do
        coefficients[i+1] := lower[i+1];
        if degree-i eq i then
            assert lower[i+1] eq sign*lower[i+1];
        end if;
        coefficients[degree-i+1] := sign*p^(degree-2*i)*lower[i+1];
    end for;
    return rec<LFUN_FORMAT | Polynomial := PolynomialRing(Integers())!coefficients,
        Degree := degree, Sign := sign, Lower := lower[2..#lower], TraceSums := sums,
        Height := height, TwistHeight := twist_height,
        TamagawaExponent := tamagawa_exponent>;
end intrinsic;
