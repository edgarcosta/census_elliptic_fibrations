import "../src/models.m": valid_model;

INVARIANT_ROW_FORMAT := recformat<Version, Q, Key, Kind, Minimality, Status,
    Degree, Sign, LowerL, Fibres, Torsion, ShaReg, AdmissibleSha, SplitMethod,
    D2, Sha2Dimensions, RefinementStatus, Message>;
ROW_HEADER := "schema_version,q,key,kind,minimality,status,degree,sign,lower_L,fibres,torsion,sha_reg,admissible_sha,split_method,d2,sha2_dimensions,refinement_status,message";

function csv_quote(s)
    error if &or[StringToCode(c) gt 127 or StringToCode(c) eq 0 : c in Eltseq(s)],
        "record cells must be ASCII without NUL";
    if &or[c in s : c in [",", "\"", "\r", "\n"]] then
        return "\"" cat Join(Split(s, "\"" : IncludeEmpty := true), "\"\"") cat "\"";
    end if;
    return s;
end function;

function csv_cells(s)
    cells := [Strings() |];
    cell := "";
    quoted := false;
    closed := false;
    i := 1;
    while i le #s do
        c := s[i];
        error if StringToCode(c) gt 127 or StringToCode(c) eq 0,
            "record cells must be ASCII without NUL";
        if quoted then
            if c eq "\"" then
                if i lt #s and s[i+1] eq "\"" then
                    cell cat:= "\"";
                    i +:= 1;
                else
                    quoted := false;
                    closed := true;
                end if;
            else
                cell cat:= c;
            end if;
        elif c eq "," then
            Append(~cells, cell);
            cell := "";
            closed := false;
        elif c eq "\"" then
            error if cell ne "" or closed, "invalid CSV quote";
            quoted := true;
        else
            error if closed or c in "\r\n", "invalid CSV cell";
            cell cat:= c;
        end if;
        i +:= 1;
    end while;
    error if quoted, "unterminated CSV quote";
    Append(~cells, cell);
    return cells;
end function;

function integer_cell(s, label, lo, hi)
    error if s eq "", label cat " must be an integer";
    negative := s[1] eq "-";
    digits := negative select s[2..#s] else s;
    error if digits eq "" or not &and[c in "0123456789" : c in Eltseq(digits)],
        label cat " must be an integer";
    n := StringToInteger(s);
    error if IntegerToString(n) ne s or n lt lo or (hi ge lo and n gt hi),
        label cat " is outside its range or is not canonical";
    return n;
end function;

function integer_list(s, label, lo, hi)
    if s eq "" then return ""; end if;
    if s eq "[]" then return [Integers() |]; end if;
    parts := Split(s, " " : IncludeEmpty := true);
    return [integer_cell(c, label, lo, hi) : c in parts];
end function;

function wire_list(values)
    if values cmpeq "" then return ""; end if;
    return #values eq 0 select "[]" else Join([IntegerToString(n) : n in values], " ");
end function;

function wire_scalar(value)
    if value cmpeq "" then return ""; end if;
    return IntegerToString(value);
end function;

function wire_rational(value)
    if value cmpeq "" then return ""; end if;
    value := Rationals()!value;
    return IntegerToString(Numerator(value)) cat "/" cat IntegerToString(Denominator(value));
end function;

function rational_cell(s)
    if s eq "" then return ""; end if;
    parts := Split(s, "/" : IncludeEmpty := true);
    error if #parts ne 2, "sha_reg must be a reduced rational";
    a := integer_cell(parts[1], "sha_reg numerator", 1, 0);
    b := integer_cell(parts[2], "sha_reg denominator", 1, 0);
    error if GCD(a, b) ne 1, "sha_reg must be a reduced rational";
    return a/b;
end function;

function fibre_cell(s)
    if s eq "" then return ""; end if;
    rows := [* *];
    if s eq "[]" then return rows; end if;
    for entry in Split(s, ";" : IncludeEmpty := true) do
        c := Split(entry, ":" : IncludeEmpty := true);
        error if #c ne 8, "fibre tuple must have eight entries";
        place := c[1] eq "inf" select "inf" else integer_cell(c[1], "place", 0, -1);
        degree := integer_cell(c[2], "place degree", 1, 255);
        error if place cmpeq "inf" and degree ne 1, "infinity must have degree one";
        error if c[3] notin ["I", "II", "III", "IV", "I*", "IV*", "III*", "II*"],
            "unknown Kodaira type";
        n := integer_cell(c[4], "fibre index", 0, 65535);
        delta := integer_cell(c[5], "discriminant valuation", 0, 65535);
        conductor := integer_cell(c[6], "conductor", 0, 65535);
        tamagawa := integer_cell(c[7], "Tamagawa number", 1, 65535);
        trace := integer_cell(c[8], "trace", -1, 1);
        Append(~rows, <place, degree, c[3], n, delta, conductor, tamagawa, trace>);
    end for;
    return rows;
end function;

function wire_fibres(rows)
    if rows cmpeq "" then return ""; end if;
    if #rows eq 0 then return "[]"; end if;
    entries := [];
    for r in rows do
        error if #r ne 8, "fibre tuple must have eight entries";
        place := r[1] cmpeq "inf" select "inf" else IntegerToString(r[1]);
        Append(~entries, Join([place, IntegerToString(r[2]), r[3]]
            cat [IntegerToString(r[i]) : i in [4..8]], ":"));
    end for;
    return Join(entries, ";");
end function;

intrinsic FibrationInvariantRow(key::RngIntElt, q::RngIntElt,
    kind::MonStgElt, minimality::MonStgElt) -> Rec
{A stage 1 invariant row with unknown arithmetic for key over GF(q), with
its generic-fibre kind and minimality. Optional unknown fields are the empty
string. See also FibrationInvariantCSV}
    require key ge 0 and key lt 2^64: "key must fit uint64";
    require q in [2, 3]: "q must be 2 or 3";
    require kind in ["smooth", "regular cusp", "nonregular"]: "unknown kind";
    require minimality in ["everywhere minimal", "nonminimal somewhere"]:
        "unknown minimality";
    supported := kind eq "smooth" and minimality eq "everywhere minimal";
    return rec<INVARIANT_ROW_FORMAT | Version := 1, Q := q, Key := key,
        Kind := kind, Minimality := minimality,
        Status := supported select "partial" else "unsupported",
        Degree := "", Sign := "", LowerL := "", Fibres := "", Torsion := "",
        ShaReg := "", AdmissibleSha := "", SplitMethod := "", D2 := "",
        Sha2Dimensions := "", RefinementStatus := supported select "" else "unsupported",
        Message := "">;
end intrinsic;

intrinsic FibrationInvariantRow(text::MonStgElt) -> Rec
{The stage 1 invariant row in one versioned CSV record text, without a trailing
record separator. This reference codec does not read release files or manifests.
See also FibrationInvariantCSV}
    c := csv_cells(text);
    require #c eq 18: "invariant CSV must have 18 columns";
    require c[1] eq "1": "unsupported schema version";
    row := FibrationInvariantRow(integer_cell(c[3], "key", 0, 2^64-1),
        integer_cell(c[2], "q", 2, 3), c[4], c[5]);
    require c[6] in ["ok", "partial", "unsupported", "error"]: "unknown status";
    row`Status := c[6];
    row`Degree := c[7] eq "" select "" else integer_cell(c[7], "degree", 0, 20);
    row`Sign := c[8] eq "" select "" else integer_cell(c[8], "sign", -1, 1);
    require row`Sign cmpeq "" or row`Sign in [-1, 1]: "sign must be -1 or 1";
    // Lower-half coefficients are unbounded signed integers.
    if c[9] eq "" or c[9] eq "[]" then
        row`LowerL := c[9] eq "" select "" else [Integers() |];
    else
        parts := Split(c[9], " " : IncludeEmpty := true);
        row`LowerL := [StringToInteger(s) : s in parts];
        require wire_list(row`LowerL) eq c[9]: "lower_L must contain canonical integers";
    end if;
    row`Fibres := fibre_cell(c[10]);
    row`Torsion := integer_list(c[11], "torsion", 2, 1);
    row`ShaReg := rational_cell(c[12]);
    row`AdmissibleSha := integer_list(c[13], "admissible_sha", 1, 0);
    require c[14] in ["", "constraints", "unresolved"]: "unknown stage 1 split method";
    row`SplitMethod := c[14];
    require c[15] eq "" and c[16] in ["", "[]"]
        and c[17] in ["", "not_needed", "unsupported", "error"]:
        "stage 1 refinement must not claim uncomputed exact or bound data";
    row`Sha2Dimensions := c[16] eq "" select "" else [Integers() |];
    row`RefinementStatus := c[17];
    row`Message := c[18];
    if not row`LowerL cmpeq "" then
        require not row`Degree cmpeq "" and #row`LowerL eq row`Degree div 2:
            "lower_L length must equal floor(degree/2)";
    end if;
    if not row`Torsion cmpeq "" then
        require &and[row`Torsion[i+1] mod row`Torsion[i] eq 0
            : i in [1..#row`Torsion-1]]: "torsion must be in divisibility order";
    end if;
    if not row`AdmissibleSha cmpeq "" then
        sha := row`AdmissibleSha;
        require sha eq Sort(Setseq(Seqset(sha))) and &and[IsSquare(n) : n in sha]:
            "admissible_sha must be sorted distinct positive squares";
    end if;
    if row`SplitMethod eq "constraints" or row`RefinementStatus eq "not_needed" then
        require not row`AdmissibleSha cmpeq "" and #row`AdmissibleSha eq 1:
            "a constraints split or not_needed refinement requires a singleton";
    end if;
    supported := row`Kind eq "smooth" and row`Minimality eq "everywhere minimal";
    if row`Status eq "ok" then
        require supported and &and[c[i] ne "" : i in [7..14]]
            and #row`AdmissibleSha gt 0: "ok requires complete stage 1 arithmetic";
    end if;
    if not supported then
        require row`Status in ["unsupported", "error"]
            and &and[c[i] eq "" : i in [7..16]]
            and row`RefinementStatus eq "unsupported":
            "other census cells have unsupported arithmetic";
    end if;
    return row;
end intrinsic;

intrinsic FibrationInvariantCSV() -> MonStgElt
{The column header of the version 1 stage 1 invariant CSV format}
    return ROW_HEADER;
end intrinsic;

intrinsic FibrationInvariantCSV(row::Rec) -> MonStgElt
{The versioned ASCII CSV encoding of a stage 1 invariant row, without a trailing
record separator. Rationals are reduced and optional unknowns remain empty cells.
See also FibrationInvariantRow}
    cells := [IntegerToString(row`Version), IntegerToString(row`Q),
        IntegerToString(row`Key), row`Kind, row`Minimality, row`Status,
        wire_scalar(row`Degree), wire_scalar(row`Sign), wire_list(row`LowerL),
        wire_fibres(row`Fibres), wire_list(row`Torsion), wire_rational(row`ShaReg),
        wire_list(row`AdmissibleSha), row`SplitMethod, wire_scalar(row`D2),
        wire_list(row`Sha2Dimensions), row`RefinementStatus, row`Message];
    text := Join([csv_quote(c) : c in cells], ",");
    checked := FibrationInvariantRow(text);
    return text;
end intrinsic;

intrinsic FibrationInvariantRow(key::RngIntElt, W::SeqEnum) -> Rec
{The stage 1 invariant row for a homogeneous census model W with key.
Only smooth everywhere minimal models receive K3 arithmetic.
See also FibrationProperties, FibrationCandidates}
    require valid_model(W): "W must have five homogeneous coefficients of degrees 2,4,6,8,12 over one supported ring";
    require key eq FibrationKey(W): "key must equal FibrationKey(W)";
    kind, minimality := FibrationProperties(W);
    q := #BaseRing(Universe(W));
    row := FibrationInvariantRow(key, q, kind, minimality);
    if row`Status eq "unsupported" then return row; end if;
    F<t> := FunctionField(GF(q));
    E := EllipticCurve([F | Evaluate(a, [t, F!1]) : a in W]);
    L := FibrationLPolynomial(E);
    fibres := FibrationLocalData(E);
    torsion := FibrationTorsion(E);
    data := FibrationInvariants(E, L`Polynomial, fibres, torsion);
    sha, regulators := FibrationCandidates(data, fibres);
    row`Degree := L`Degree;
    row`Sign := L`Sign;
    row`LowerL := L`Lower;
    row`Fibres := [* <r`Place, r`Degree, r`Type, r`Index,
        r`DiscriminantValuation, r`Conductor, r`Tamagawa, r`Trace> : r in fibres *];
    row`Torsion := torsion;
    row`ShaReg := data`ShaRegulatorProduct;
    row`AdmissibleSha := sha;
    row`SplitMethod := #sha eq 1 select "constraints" else "unresolved";
    row`RefinementStatus := #sha eq 1 select "not_needed"
        else q eq 2 select "unsupported" else "";
    row`Status := #sha gt 0 select "ok" else "error";
    if #sha eq 0 then row`Message := "no admissible Sha order"; end if;
    return row;
end intrinsic;
