function published_number(text, document, row)
    words := ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten"];
    if text in words then return Index(words, text)-1; end if;
    found, match, powers := Regexp("^([0-9]+)\\^([0-9]+)$", text);
    if found then return StringToInteger(powers[1])^StringToInteger(powers[2]); end if;
    digits := Join(Split(text, ","), "");
    error if #digits eq 0 or not &and[c in "0123456789" : c in Eltseq(digits)],
        Sprintf("MISMATCH: %o, %o: cannot parse count %o", document, row, text);
    return StringToInteger(digits);
end function;

function table_columns(line, document, row)
    line := Join(Split(line, "\t\r"), " ");
    error if #line lt 2 or line[1] ne "|" or line[#line] ne "|",
        Sprintf("MISMATCH: %o, %o: cannot parse table row", document, row);
    return [Join(Split(c, " "), " ") : c in Split(line[2..#line-1], "|")];
end function;

procedure check_published_table(document, header, labels, counts)
    lines := Split(Read(document), "\n");
    indices := [i : i in [1..#lines] |
        Join(Split(lines[i], " \t\r"), " ") eq header];
    error if #indices ne 1, Sprintf(
        "MISMATCH: %o: expected one census table headed %o, found %o",
        document, header, #indices);
    first := indices[1]+2;
    error if first+#labels-1 gt #lines,
        Sprintf("MISMATCH: %o: incomplete census table", document);
    separator := table_columns(lines[first-1], document, "header separator");
    error if #separator ne #labels[1]+2 or not &and[
        "-" in c and &and[x in "-:" : x in Eltseq(c)] : c in separator],
        Sprintf("MISMATCH: %o: cannot parse census table separator", document);
    for i in [1..#labels] do
        row := Join(labels[i], ", ");
        columns := table_columns(lines[first+i-1], document, row);
        error if #columns ne #labels[i]+2, Sprintf(
            "MISMATCH: %o, %o: expected %o columns, found %o",
            document, row, #labels[i]+2, #columns);
        error if columns[1..#labels[i]] ne labels[i], Sprintf(
            "MISMATCH: %o, %o: missing or unrecognized census row %o",
            document, row, columns);
        for j in [1, 2] do
            actual := published_number(columns[#labels[i]+j], document, row);
            error if actual ne counts[j][i], Sprintf(
                "MISMATCH: %o, %o, GF(%o): document %o, NumberOfFibrations %o",
                document, row, j+1, actual, counts[j][i]);
        end for;
    end for;
    last := first+#labels;
    error if last le #lines and "|" in lines[last],
        Sprintf("MISMATCH: %o: unexpected census table row %o", document, lines[last]);
end procedure;

procedure check_published_prefixes(document, p, expected)
    marker := Sprintf("prefixes over F%o.", p);
    lines := [line : line in Split(Read(document), "\n") | marker in line];
    error if #lines ne 1, Sprintf(
        "MISMATCH: %o, F%o prefixes: expected one count, found %o", document, p, #lines);
    found, match, groups := Regexp(Sprintf("([0-9,]+) prefixes over F%o[.]", p), lines[1]);
    error if not found, Sprintf(
        "MISMATCH: %o, F%o prefixes: cannot parse count", document, p);
    actual := published_number(groups[1], document, marker);
    error if actual ne expected, Sprintf(
        "MISMATCH: %o, F%o prefixes: document %o, FibrationPrefixes %o",
        document, p, actual, expected);
end procedure;

procedure check_published_statement(document, pattern, expected, source)
    text := Join(Split(Read(document), " \t\r\n"), " ");
    found, match, groups := Regexp(pattern, text);
    error if not found, Sprintf(
        "MISMATCH: %o: missing statement or block matching %o", document, pattern);
    error if #groups ne #expected, Sprintf(
        "MISMATCH: %o, %o: expected %o counts, found %o",
        document, match, #expected, #groups);
    for i in [1..#expected] do
        actual := published_number(groups[i], document, match);
        error if actual ne expected[i], Sprintf(
            "MISMATCH: %o, %o: count %o: document %o, %o %o",
            document, match, i, actual, source, expected[i]);
    end for;
end procedure;

// Pins census tables, example outputs, listing size, prefix and equation counts
// in the README, mathematics and interface docs. Missing statements must fail.
procedure test_published_counts()
    labels := [];
    mathematics_labels := [];
    kinds := ["smooth", "regular cusp", "nonregular"];
    titles := ["Smooth", "Regular cusp", "Nonregular"];
    for i in [1..#kinds] do
        for minimality in ["everywhere minimal", "nonminimal somewhere"] do
            Append(~labels, [kinds[i], minimality]);
            Append(~mathematics_labels, [titles[i] cat ", " cat minimality]);
        end for;
    end for;
    counts := [];
    prefix_counts := [];
    for p in [2, 3] do
        R<s, t> := PolynomialRing(GF(p), 2);
        n, cells := NumberOfFibrations(R);
        Append(~counts, [cells[<label[1], label[2]>] : label in labels] cat [n]);
        Append(~prefix_counts, #FibrationPrefixes(R));
        check_published_prefixes("docs/mathematics.md", p, prefix_counts[#prefix_counts]);
        check_published_statement("docs/interface.md",
            Sprintf("GF[(]%o[)]: ([0-9,]+) prefixes, ([0-9,]+) quasi-elliptic K3 orbits", p),
            [prefix_counts[p-1], cells[<"regular cusp", "everywhere minimal">]],
            "FibrationPrefixes and NumberOfFibrations");
        // equation_counts.m checks these totals through NumberOfEquations.
        check_published_statement("docs/interface.md",
            Sprintf("`(%o\\^[0-9]+)` equations with no selections[.]", p),
            [p eq 2 select 2^32 else 3^27], "normalized equation total");
    end for;
    check_published_statement("docs/interface.md",
        "The tables have ([0-9,]+) and ([0-9,]+) entries, respectively[.]",
        prefix_counts, "FibrationPrefixes");
    check_published_table("README.md",
        "| Generic-fibre kind | Minimality | GF(2) | GF(3) |", labels cat [["All", ""]], counts);
    check_published_table("docs/mathematics.md", "| Quantity | F2 | F3 |",
        [["All orbits"]] cat mathematics_labels, [[c[7]] cat c[1..6] : c in counts]);
    check_published_statement("README.md",
        "```text ([0-9,]+) ([0-9,]+) ([0-9,]+) ```",
        [counts[2][7], counts[2][1], counts[2][3]], "NumberOfFibrations");

    R<s, t> := PolynomialRing(GF(3), 2);
    prefixes := FibrationPrefixes(R);
    i := Index(prefixes, <R!0, R!0>);
    keys := FibrationRepresentativeKeys(R, i : GenericFibreKinds := ["regular cusp"],
        Minimalities := ["everywhere minimal"]);
    check_published_statement("README.md",
        "```text ([0-9,]+) regular cusp everywhere minimal ```",
        [#keys], "FibrationRepresentativeKeys");
    // The binary format in docs/keys.md stores eight bytes per key.
    check_published_statement("README.md",
        "that is ([0-9,]+) keys, occupying ([0-9,]+) bytes at eight bytes per key, excluding metadata[.]",
        [counts[2][7], 8*counts[2][7]], "NumberOfFibrations and bytes per key");
end procedure;

test_published_counts();
