import "tests/support/_support.m": expect_error, shell_quote;

function curve_from_key(p, key)
    R<s, u> := PolynomialRing(GF(p), 2);
    F<t> := FunctionField(GF(p));
    return EllipticCurve([F | Evaluate(a, [t, F!1])
        : a in aInvariantsFromKey(R, key)]);
end function;

function independent_reference(E)
    p := Characteristic(BaseRing(E));
    inputpath, inputfile := TemporaryFile();
    outputpath, outputfile := TemporaryFile();
    logpath, logfile := TemporaryFile();
    delete inputfile;
    delete outputfile;
    delete logfile;
    a := [SequenceToInteger([Integers()!c : c in Coefficients(Numerator(b))], p)
        : b in aInvariants(E)];
    Write(inputpath, Join([Sprint(c) : c in a], " "));
    try
        command := "(ulimit -t 120; exec \"${MAGMA:-magma}\" -bn p:=" cat Sprint(p)
            cat " inputpath:=" cat shell_quote(inputpath)
            cat " outputpath:=" cat shell_quote(outputpath)
            cat " tests/support/magma/lfun_reference.m) > " cat shell_quote(logpath) cat " 2>&1";
        status := System(command);
        error if status div 256 in [137, 152] or status mod 128 in [9, 24],
            "required independent LFunction reference timed out (120 CPU seconds)";
        error if status ne 0 or #Read(outputpath) eq 0, Read(logpath);
        L := PolynomialRing(Integers())![StringToInteger(c) : c in Split(Read(outputpath), " \n")];
    catch e
        Remove(inputpath);
        Remove(outputpath);
        Remove(logpath);
        error e;
    end try;
    Remove(inputpath);
    Remove(outputpath);
    Remove(logpath);
    return L;
end function;

// Built-in LFunction shares Magma's local-fibre routines with the package.
// Direct fibre counts and hand reductions certify the known examples.
procedure check_polynomial(E, expected_sign, rank_zero, open_sign)
    result := FibrationLPolynomial(E);
    reference := independent_reference(E);
    p := Characteristic(BaseRing(E));
    assert result`Polynomial eq reference;
    assert result`Degree eq Degree(reference);
    assert result`Sign eq expected_sign;
    assert LeadingCoefficient(reference) eq expected_sign*p^Degree(reference);
    assert (Evaluate(reference, 1/p) ne 0) eq rank_zero;
    D := result`Degree;
    c := Coefficients(result`Polynomial);
    assert #result`Lower eq D div 2;
    assert result`Lower eq c[2..1+D div 2];
    for i in [0..D div 2] do
        assert c[D-i+1] eq expected_sign*p^(D-2*i)*c[i+1];
    end for;
    if open_sign then
        assert &and[Coefficient(reference, i) eq 0 : i in [D div 2..12]];
    end if;
    // Newton sums from the built-in cross-check pin every counted extension.
    for n in [1..#result`TraceSums] do
        assert n*c[n+1] eq &+[result`TraceSums[j]*c[n-j+1] : j in [1..n]];
    end for;
    F := BaseRing(E);
    t := F.1;
    weights := [2, 4, 6, 8, 12];
    a := aInvariants(E);
    reciprocal := EllipticCurve([F | t^weights[i]*Evaluate(a[i], 1/t) : i in [1..5]]);
    other := FibrationLPolynomial(reciprocal);
    assert other`Polynomial eq reference;
    assert other`Sign eq result`Sign and other`Lower eq result`Lower;
    assert other`TraceSums eq result`TraceSums;
end procedure;

// Formerly open signs include both characteristics, both signs and rank zero.
procedure test_former_open_signs()
    check_polynomial(curve_from_key(2, 4884480), 1, true, true);
    check_polynomial(curve_from_key(2, 11821184), -1, false, true);
    check_polynomial(curve_from_key(3, 372699133912), -1, false, true);
    check_polynomial(curve_from_key(3, 320009441124), 1, true, true);
end procedure;

// Wild conductors and minimized twists determine different degrees and signs.
procedure test_wild_and_ordinary_fibres()
    F<t> := FunctionField(GF(2));
    for d in [6, 7] do
        E := EllipticCurve([F | 0, 0, 1, t^d, 0]);
        check_polynomial(E, d eq 6 select 1 else -1, false, false);
        result := FibrationLPolynomial(E);
        assert <result`Height, result`TwistHeight, result`TamagawaExponent>
            eq (d eq 6 select <2, 1, -1> else <2, 4, 1>);
        assert result`Degree eq (d eq 6 select 7 else 19);
    end for;
    check_polynomial(EllipticCurve([F | t, 0, t^6+1, 0, 0]), -1, false, false);
    F<t> := FunctionField(GF(3));
    E := EllipticCurve([F | 0, 0, 0, t^8, t+1]);
    check_polynomial(E, 1, false, false);
    result := FibrationLPolynomial(E);
    assert <result`Height, result`TwistHeight, result`TamagawaExponent>
        eq <2, 2, 0>;
    assert result`Degree eq 20;
end procedure;

// Coordinate changes preserve the polynomial even with nonzero a2, a1 or a3.
procedure test_unnormalized_equations()
    F<t> := FunctionField(GF(2));
    check_polynomial(EllipticCurve([F | t, t, t^6+t^2+1, t^2, t^3]),
        -1, false, false);
    F<t> := FunctionField(GF(3));
    check_polynomial(EllipticCurve([F | 1, t^4, t, t^4+t, -t^2]),
        -1, false, false);
end procedure;

// The sign proof applies to the supported K3 population and polynomial chart.
procedure test_lfun_domain()
    F<t> := FunctionField(GF(4));
    E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
    expect_error(procedure() result := FibrationLPolynomial(E); end procedure,
        "constant field must be GF(2) or GF(3)");
    F<t> := FunctionField(GF(2));
    for a in [[F | 0, 0, 1, 0, 0], [F | 0, 0, t^6, t^8, 0]] do
        E := EllipticCurve(a);
        expect_error(procedure() result := FibrationLPolynomial(E); end procedure,
            "equation must define an everywhere minimal elliptic K3");
    end for;
    for a4 in [F | 1/t, t^9] do
        E := EllipticCurve([F | 0, 0, 1, a4, 0]);
        expect_error(procedure() result := FibrationLPolynomial(E); end procedure,
            "coefficients must be polynomials");
    end for;
end procedure;

test_former_open_signs();
test_wild_and_ordinary_fibres();
test_unnormalized_equations();
test_lfun_domain();
