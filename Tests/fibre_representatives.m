import "Tests/_support.m": expect_error;

// Pins listings over the cuspidal prefixes, including empty smooth selections.
procedure test_selected_representatives()
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        keys := FibrationRepresentativeKeys(R,1);
        assert #keys eq (p eq 2 select 720 else 315);
        assert #Seqset(keys) eq #keys;
        assert IsEmpty(FibrationRepresentativeKeys(R,1 : GenericFibreKinds := ["smooth"]));
        for kind in ["regular cusp","nonregular"] do
            for minimality in ["everywhere minimal","nonminimal somewhere"] do
                selected := FibrationRepresentativeKeys(R,1 : GenericFibreKinds := [kind], Minimalities := [minimality]);
                n,cells,rows := NumberOfFibrations(R : PrefixIndices := [1],
                    GenericFibreKinds := [kind], Minimalities := [minimality]);
                if n eq 0 then assert IsEmpty(selected); continue; end if;
                assert #selected gt 0 and #selected eq n;
                for key in selected[1..Min(#selected,3)] do
                    k,m := FibrationProperties(aInvariantsFromKey(R,key));
                    assert <k,m> eq <kind,minimality>;
                end for;
            end for;
        end for;
        smooth := FibrationRepresentativeKeys(R,2 : GenericFibreKinds := ["smooth"]);
        assert #smooth eq (p eq 2 select 4296 else 1143);
        assert #Seqset(smooth) eq #smooth;
        for minimality in ["everywhere minimal","nonminimal somewhere"] do
            selected := FibrationRepresentativeKeys(R,2 : GenericFibreKinds := ["smooth"], Minimalities := [minimality]);
            n,cells,rows := NumberOfFibrations(R : PrefixIndices := [2],
                GenericFibreKinds := ["smooth"], Minimalities := [minimality]);
            assert #selected gt 0 and #selected eq n;
            for key in selected[1..Min(#selected,3)] do
                k,m := FibrationProperties(aInvariantsFromKey(R,key));
                assert k eq "smooth" and m eq minimality;
            end for;
        end for;
        expect_error(procedure() x := FibrationRepresentativeKeys(R,0); end procedure,"valid prefix indices");
        expect_error(procedure() n,c,r := NumberOfFibrations(R : PrefixIndices := [1,1]); end procedure,"PrefixIndices must be distinct");
        expect_error(procedure() n,c,r := NumberOfFibrations(R : GenericFibreKinds := ["Smooth"]); end procedure,"GenericFibreKinds must contain distinct values");
        expect_error(procedure() n,c,r := NumberOfFibrations(R : GenericFibreKinds := ["smooth","smooth"]); end procedure,"GenericFibreKinds must contain distinct values");
        expect_error(procedure() n,c,r := NumberOfFibrations(R : Minimalities := ["minimal"]); end procedure,"Minimalities must contain distinct values");
        expect_error(procedure() n,c,r := NumberOfFibrations(R : Minimalities := ["everywhere minimal","everywhere minimal"]); end procedure,"Minimalities must contain distinct values");
    end for;
end procedure;

test_selected_representatives();
