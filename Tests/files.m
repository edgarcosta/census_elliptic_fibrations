import "Tests/_support.m": expect_error;

procedure remove_outputs(stem, indices)
    for path in [ stem cat ".csv" ] cat &cat[ [ Sprintf("%o_%o.bin",stem,IntegerToString(10000+i)[2..5]),
        Sprintf("%o_%o.stats",stem,IntegerToString(10000+i)[2..5]) ] : i in indices ] do
        if FileExists(path) then Remove(path); end if;
    end for;
end procedure;

// Pins exact metadata, zero cells and byte order against the cuspidal counts.
procedure test_representative_files()
    stem,handle := TemporaryFile();
    delete handle;
    Remove(stem);
    try
        for p in [2,3] do
            R<s,t> := PolynomialRing(GF(p),2);
            n := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1]);
            expected := p eq 2 select [0,0,680,38,0,2] else [0,0,294,20,0,1];
            assert n eq &+expected;
            keys := FibrationRepresentativeKeys(R,1);
            F := Open(stem cat "_0001.bin","rb");
            bytes := ReadBytes(F);
            delete F;
            assert #bytes eq 8*n and n gt 0;
            for i in [1..n] do
                assert &+[ bytes[8*(i-1)+j+1]*256^j : j in [0..7] ] eq keys[i];
            end for;
            degrees := p eq 2 select [2,6] else [4,8];
            text := Sprintf("file_format_version=1\nkey_encoding_version=1\ncharacteristic=%o\nprefix_index=1\nprefix_form_1_degree=%o\nprefix_form_1_key=0\nprefix_form_2_degree=%o\nprefix_form_2_key=0\ngeneric_fibre_kinds=smooth|regular cusp|nonregular\nminimalities=everywhere minimal|nonminimal somewhere\nnumber_of_keys=%o\n",p,degrees[1],degrees[2],n);
            csv := "file_format_version,key_encoding_version,characteristic,prefix_index,prefix_form_1_degree,prefix_form_1_key,prefix_form_2_degree,prefix_form_2_key,generic_fibre_kind,minimality,number_of_keys\n";
            j := 0;
            for kind in ["smooth","regular cusp","nonregular"] do
                for minimality in ["everywhere minimal","nonminimal somewhere"] do
                    j +:= 1;
                    text cat:= Sprintf("cell=%o|%o|%o\n",kind,minimality,expected[j]);
                    csv cat:= Sprintf("1,1,%o,1,%o,0,%o,0,%o,%o,%o\n",p,degrees[1],degrees[2],kind,minimality,expected[j]);
                end for;
            end for;
            assert Read(stem cat "_0001.stats") eq text cat "complete=1\n";
            assert Read(stem cat ".csv") eq csv;
            expect_error(procedure() x := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1]); end procedure,"An output path already exists");
            remove_outputs(stem,[1]);
            empty := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1], GenericFibreKinds := ["smooth"]);
            assert empty eq 0;
            F := Open(stem cat "_0001.bin","rb");
            assert IsEmpty(ReadBytes(F));
            delete F;
            assert "number_of_keys=0\n" in Read(stem cat "_0001.stats");
            remove_outputs(stem,[1]);
        end for;
        R<s,t> := PolynomialRing(GF(3),2);
        n := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1],GenericFibreKinds := ["nonregular"]);
        assert n eq 1;
        F := Open(stem cat "_0001.bin","rb");
        assert ReadBytes(F) eq [0,0,0,0,0,0,0,0];
        delete F;
        remove_outputs(stem,[1]);
        n := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [2,1],
            GenericFibreKinds := ["regular cusp"],
            Minimalities := ["nonminimal somewhere","everywhere minimal"]);
        assert n eq 314;
        lines := Split(Read(stem cat ".csv"),"\n");
        assert #lines eq 5;
        assert "1,1,3,2," eq lines[2][1..8] and "1,1,3,1," eq lines[4][1..8];
        assert "minimalities=everywhere minimal|nonminimal somewhere\n" in Read(stem cat "_0001.stats");
        remove_outputs(stem,[1,2]);
        expect_error(procedure() x := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1],Minimalities := ["unknown"]); end procedure,"Minimalities must contain distinct values");
        assert not FileExists(stem cat ".csv");
        Write(stem cat "_0002.stats","incomplete");
        expect_error(procedure() x := WriteFibrationRepresentativeKeys(R,stem : PrefixIndices := [1,2]); end procedure,"An output path already exists");
        assert not FileExists(stem cat ".csv") and not FileExists(stem cat "_0001.bin");
    catch e
        remove_outputs(stem,[1,2]);
        error e;
    end try;
    remove_outputs(stem,[1,2]);
end procedure;

test_representative_files();
