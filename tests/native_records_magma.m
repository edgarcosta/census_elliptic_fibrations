import "tests/support/_support.m": expect_error;
import "tests/support/_native_model_samples.m": model_samples;

// Literal CSV pins field order, canonical rationals and the two empty meanings.
procedure test_exact_record_codec()
    row := FibrationInvariantRow(2^64-1, 3, "smooth", "everywhere minimal");
    row`Degree := 4;
    row`Sign := -1;
    row`LowerL := [-(2^80+1), 0];
    row`Fibres := [* <81, 4, "I", 6, 6, 1, 2, -1>,
        <"inf", 1, "II*", 0, 12, 4, 1, 0> *];
    row`Torsion := [2, 4];
    row`ShaReg := (6*(2^90+1))/8;
    row`AdmissibleSha := [1, 4, 9^20];
    row`SplitMethod := "unresolved";
    row`Sha2Dimensions := [];
    row`Message := "trace, \"pending\"\nretry";
    expected := "1,3,18446744073709551615,smooth,everywhere minimal,partial,4,-1,"
        cat "-1208925819614629174706177 0,81:4:I:6:6:1:2:-1;inf:1:II*:0:12:4:1:0,"
        cat "2 4,3713820117856140824697372675/4,1 4 12157665459056928801,unresolved,,[],,"
        cat "\"trace, \"\"pending\"\"\nretry\"";
    text := FibrationInvariantCSV(row);
    assert text eq expected;
    copy := FibrationInvariantRow(text);
    assert copy`Key eq 2^64-1 and copy`LowerL eq row`LowerL;
    assert copy`Fibres eq row`Fibres and copy`Torsion eq [2, 4];
    assert copy`ShaReg eq row`ShaReg and copy`AdmissibleSha eq row`AdmissibleSha;
    assert copy`D2 cmpeq "" and copy`Sha2Dimensions eq [];
    assert copy`RefinementStatus eq "" and copy`Message eq row`Message;
    assert FibrationInvariantCSV(copy) eq text;
    assert FibrationInvariantCSV() eq "schema_version,q,key,kind,minimality,status,degree,sign,lower_L,fibres,torsion,sha_reg,admissible_sha,split_method,d2,sha2_dimensions,refinement_status,message";
    row`LowerL := [];
    row`Degree := 0;
    row`Sign := 1;
    row`Fibres := [* *];
    row`Torsion := [];
    row`ShaReg := 3;
    row`AdmissibleSha := [];
    row`Status := "error";
    row`RefinementStatus := "error";
    copy := FibrationInvariantRow(FibrationInvariantCSV(row));
    assert copy`LowerL eq [] and #copy`Fibres eq 0 and copy`Torsion eq [];
    assert copy`ShaReg eq 3 and copy`AdmissibleSha eq [];
    assert copy`Status eq "error" and copy`RefinementStatus eq "error";
end procedure;

// Classification covers every cell without applying K3 arithmetic to the other five.
procedure test_six_cell_records()
    for p in [2, 3] do
        R := PolynomialRing(GF(p), 2);
        samples := model_samples(R);
        for i in [1..6] do
            W := samples[i];
            key := FibrationKey(W);
            row := FibrationInvariantRow(key, W);
            kind, minimality := FibrationProperties(W);
            assert row`Kind eq kind and row`Minimality eq minimality;
            copy := FibrationInvariantRow(FibrationInvariantCSV(row));
            assert copy`Kind eq kind and copy`Minimality eq minimality;
            if i eq 1 then
                assert copy`Status eq "ok" and copy`ShaReg gt 0;
            else
                assert copy`Status eq "unsupported";
                assert copy`Degree cmpeq "" and copy`Torsion cmpeq "";
                assert copy`AdmissibleSha cmpeq "";
                assert copy`RefinementStatus eq "unsupported";
            end if;
        end for;
    end for;
end procedure;

// Fixed curves pin the compressed L and ambiguous lists against packaged arithmetic.
procedure test_computed_records()
    for p in [2, 3] do
        R<s, u> := PolynomialRing(GF(p), 2);
        F<t> := FunctionField(GF(p));
        if p eq 2 then
            W := [R | 0, 0, u^6, s^7*u, 0];
            E := EllipticCurve([F | 0, 0, 1, t^7, 0]);
            sha := [1, 4];
            product := 14;
        else
            W := [R | 0, 0, 0, u^8, s^7*u^5];
            E := EllipticCurve([F | 0, 0, 0, 1, t^7]);
            sha := [1, 4, 9];
            product := 108;
        end if;
        row := FibrationInvariantRow(FibrationKey(W), W);
        row := FibrationInvariantRow(FibrationInvariantCSV(row));
        reference := FibrationLPolynomial(E);
        assert row`Degree eq reference`Degree and row`Sign eq reference`Sign;
        assert row`LowerL eq reference`Lower;
        P<T> := PolynomialRing(Integers());
        D := row`Degree;
        L := P!([1] cat row`LowerL);
        for i in [0..(D-1) div 2] do
            L +:= row`Sign*p^(D-2*i)*Coefficient(L, i)*T^(D-i);
        end for;
        assert L eq P!LFunction(E);
        assert row`ShaReg eq product and row`AdmissibleSha eq sha;
        assert row`SplitMethod eq "unresolved";
        assert row`D2 cmpeq "" and row`Sha2Dimensions cmpeq "";
        assert row`RefinementStatus eq (p eq 2 select "unsupported" else "");
        assert #row`Fibres eq #FibrationLocalData(E);
        assert row`Torsion eq FibrationTorsion(E);
    end for;
    R<s, u> := PolynomialRing(GF(2), 2);
    W := [R | s^2, 0, u^6, 0, 0];
    row := FibrationInvariantRow(FibrationKey(W), W);
    assert row`AdmissibleSha eq [1] and row`ShaReg eq 1;
    assert row`SplitMethod eq "constraints" and row`RefinementStatus eq "not_needed";
    R<s, u> := PolynomialRing(GF(3), 2);
    W := [R | 0, s^4+u^4, 0, s^4*u^4, 0];
    row := FibrationInvariantRow(FibrationKey(W), W);
    row := FibrationInvariantRow(FibrationInvariantCSV(row));
    assert row`Torsion eq [2, 4] and row`ShaReg eq 3/4;
    assert row`AdmissibleSha eq [1] and row`SplitMethod eq "constraints";
end procedure;

// Invalid wire values and unsupported versions must not acquire exact meanings.
procedure test_record_rejections()
    R<s, u> := PolynomialRing(GF(2), 2);
    W := [R | s^2, 0, u^6, 0, 0];
    key := FibrationKey(W)+1;
    expect_error(procedure() row := FibrationInvariantRow(key, W); end procedure,
        "key must equal FibrationKey(W)");
    row := FibrationInvariantRow(0, 3, "smooth", "everywhere minimal");
    row`RefinementStatus := "exact";
    expect_error(procedure() text := FibrationInvariantCSV(row); end procedure,
        "stage 1 refinement");
    row`RefinementStatus := "";
    row`Key := 2^64;
    expect_error(procedure() text := FibrationInvariantCSV(row); end procedure,
        "key");
    row`Key := 0;
    row`Message := CodeToString(233);
    expect_error(procedure() text := FibrationInvariantCSV(row); end procedure,
        "ASCII");
    cells := ["2", "3", "0", "smooth", "everywhere minimal", "partial"]
        cat ["" : i in [1..12]];
    expect_error(procedure() copy := FibrationInvariantRow(Join(cells, ",")); end procedure,
        "schema version");
    cells[1] := "1";
    cells[12] := "6/8";
    expect_error(procedure() copy := FibrationInvariantRow(Join(cells, ",")); end procedure,
        "reduced rational");
    cells[12] := "1/1";
    cells[3] := "01";
    expect_error(procedure() copy := FibrationInvariantRow(Join(cells, ",")); end procedure,
        "key");
    for claim in [<2, "partial", [1], "unresolved", "">,
        <3, "partial", [1], "constraints", "unsupported">,
        <2, "ok", [1], "constraints", "">,
        <3, "partial", [1, 4], "unresolved", "unsupported">,
        <2, "ok", [1, 4], "unresolved", "">] do
        q, status, sha, split, refinement := Explode(claim);
        row := FibrationInvariantRow(0, q, "smooth", "everywhere minimal");
        row`Status := status;
        row`Degree := 0;
        row`Sign := 1;
        row`LowerL := [];
        row`Fibres := [**];
        row`Torsion := [];
        row`ShaReg := 1;
        row`AdmissibleSha := sha;
        row`SplitMethod := split;
        row`RefinementStatus := refinement;
        cells := ["1", IntegerToString(q), "0", "smooth", "everywhere minimal",
            status, "0", "1", "[]", "[]", "[]", "1/1",
            Join([IntegerToString(n) : n in sha], " "), split, "", "", refinement, ""];
        expect_error(procedure() copy := FibrationInvariantRow(Join(cells, ",")); end procedure,
            "candidate");
        expect_error(procedure() text := FibrationInvariantCSV(row); end procedure,
            "candidate");
    end for;
    for q in [2, 3], sha in [[1], [1, 4]] do
        row := FibrationInvariantRow(0, q, "smooth", "everywhere minimal");
        row`AdmissibleSha := sha;
        copy := FibrationInvariantRow(FibrationInvariantCSV(row));
        assert copy`SplitMethod eq "" and copy`RefinementStatus eq "";
        row`Status := "error";
        row`SplitMethod := #sha eq 1 select "constraints" else "unresolved";
        row`RefinementStatus := "error";
        copy := FibrationInvariantRow(FibrationInvariantCSV(row));
        assert copy`Status eq "error" and copy`RefinementStatus eq "error";
    end for;
    expect_error(procedure() copy := FibrationInvariantRow("2,3,0"); end procedure,
        "columns");
    expect_error(procedure() copy := FibrationInvariantRow("\"unterminated"); end procedure,
        "CSV");
end procedure;

test_exact_record_codec();
test_six_cell_records();
test_computed_records();
test_record_rejections();
