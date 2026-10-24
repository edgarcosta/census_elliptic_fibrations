// Pins the direct square and cube criteria, separately from minimality.
procedure test_generic_fibre_properties()
    for p in [2,3] do
        R<s,t> := PolynomialRing(GF(p),2);
        smooth := p eq 2 select [ R | 0,0,t^6,0,s^12 ] else [ R | 0,0,0,t^8,s^12 ];
        examples := [ <smooth,"smooth">,
            <[ R | 0,0,0,0,s*t^11 ],"regular cusp">,
            <[ R | 0,0,0,0,s^12 ],"nonregular"> ];
        nodal := p eq 2 select [ R | s*t,s^4+t^4,0,0,0 ]
            else [ R | 0,s*t*(s^2+t^2),0,0,0 ];
        Append(~examples,<nodal,"nonregular">);
        if p eq 2 then Append(~examples,<[ R | 0,0,0,s^7*t,0 ],"regular cusp">); end if;
        assert #examples gt 0;
        for pair in examples do
            kind,minimality := FibrationProperties(pair[1]);
            assert kind eq pair[2];
            ok,certificate := IsEverywhereMinimal(pair[1]);
            assert (minimality eq "everywhere minimal") eq ok;
        end for;
    end for;
end procedure;

test_generic_fibre_properties();
