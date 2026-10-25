function published_number(text, document, row)
    digits := Join(Split(text,","),"");
    error if #digits eq 0 or not &and[ c in "0123456789" : c in Eltseq(digits) ],
        Sprintf("MISMATCH: %o, %o: cannot parse count %o",document,row,text);
    return StringToInteger(digits);
end function;

function table_columns(line, document, row)
    line := Join(Split(line,"\t\r")," ");
    error if #line lt 2 or line[1] ne "|" or line[#line] ne "|",
        Sprintf("MISMATCH: %o, %o: cannot parse table row",document,row);
    return [ Join(Split(c," ")," ") : c in Split(line[2..#line-1],"|") ];
end function;

procedure check_published_table(document, header, labels, counts)
    lines := Split(Read(document),"\n");
    indices := [ i : i in [1..#lines] |
        Join(Split(lines[i]," \t\r")," ") eq header ];
    error if #indices ne 1, Sprintf(
        "MISMATCH: %o: expected one census table headed %o, found %o",
        document,header,#indices);
    first := indices[1]+2;
    error if first+#labels-1 gt #lines,
        Sprintf("MISMATCH: %o: incomplete census table",document);
    separator := table_columns(lines[first-1],document,"header separator");
    error if #separator ne #labels[1]+2 or not &and[
        "-" in c and &and[ x in "-:" : x in Eltseq(c) ] : c in separator ],
        Sprintf("MISMATCH: %o: cannot parse census table separator",document);
    for i in [1..#labels] do
        row := Join(labels[i],", ");
        columns := table_columns(lines[first+i-1],document,row);
        error if #columns ne #labels[i]+2, Sprintf(
            "MISMATCH: %o, %o: expected %o columns, found %o",
            document,row,#labels[i]+2,#columns);
        error if columns[1..#labels[i]] ne labels[i], Sprintf(
            "MISMATCH: %o, %o: missing or unrecognized census row %o",
            document,row,columns);
        for j in [1,2] do
            actual := published_number(columns[#labels[i]+j],document,row);
            error if actual ne counts[j][i], Sprintf(
                "MISMATCH: %o, %o, GF(%o): document %o, NumberOfFibrations %o",
                document,row,j+1,actual,counts[j][i]);
        end for;
    end for;
    last := first+#labels;
    error if last le #lines and "|" in lines[last],
        Sprintf("MISMATCH: %o: unexpected census table row %o",document,lines[last]);
end procedure;

procedure check_published_prefixes(document, p, expected)
    marker := Sprintf("prefixes over F%o.",p);
    lines := [ line : line in Split(Read(document),"\n") | marker in line ];
    error if #lines ne 1, Sprintf(
        "MISMATCH: %o, F%o prefixes: expected one count, found %o",document,p,#lines);
    found,match,groups := Regexp(Sprintf("([0-9,]+) prefixes over F%o[.]",p),lines[1]);
    error if not found, Sprintf(
        "MISMATCH: %o, F%o prefixes: cannot parse count",document,p);
    actual := published_number(groups[1],document,marker);
    error if actual ne expected, Sprintf(
        "MISMATCH: %o, F%o prefixes: document %o, FibrationPrefixes %o",
        document,p,actual,expected);
end procedure;

// Pins both published census tables and the stated prefix counts to the public API.
// Requiring their headers and every row prevents layout edits from skipping checks.
procedure test_published_counts()
    labels := [];
    mathematics_labels := [];
    kinds := ["smooth","regular cusp","nonregular"];
    titles := ["Smooth","Regular cusp","Nonregular"];
    for i in [1..#kinds] do
        for minimality in ["everywhere minimal","nonminimal somewhere"] do
            Append(~labels,[kinds[i],minimality]);
            Append(~mathematics_labels,[titles[i] cat ", " cat minimality]);
        end for;
    end for;
    counts := [];
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        n,cells := NumberOfFibrations(R);
        Append(~counts,[ cells[<label[1],label[2]>] : label in labels ] cat [n]);
        check_published_prefixes("docs/mathematics.md",p,#FibrationPrefixes(R));
    end for;
    check_published_table("README.md",
        "| Generic-fibre kind | Minimality | GF(2) | GF(3) |",labels cat [["All",""]],counts);
    check_published_table("docs/mathematics.md","| Quantity | F2 | F3 |",
        [["All orbits"]] cat mathematics_labels,[ [c[7]] cat c[1..6] : c in counts ]);
end procedure;

test_published_counts();
