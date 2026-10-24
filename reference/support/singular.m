// Moving the unique node to the origin leaves only small spaces of forms.
function reference_singular_models(R)
    s := R.1; t := R.2;
    if Characteristic(R) eq 2 then
        keys := FibrationFiberOrbitKeysF2(<R!0,R!0> : GenericSmoothOnly := false);
        models := [ FibrationFromKeyF2(R,key) : key in keys ];
        cusp_count := #models;
        nodekeys := {};
        for a1 in [s^2,s*t,s^2+s*t+t^2], n in [0..31] do
            _,key := CanonicalFibrationF2([ R | a1,EF2FormFromMask(R,4,n),0,0,0 ]);
            Include(~nodekeys,key);
        end for;
        models cat:= [ FibrationFromKeyF2(R,key) : key in nodekeys ];
    else
        P := <R!0,R!0,[* <M,R!0> : M in EF3_BaseMatrices(GF(3)) *]>;
        keys := FibrationFiberOrbitKeysF3(P : GenericSmoothOnly := false);
        triples := [ FibrationFromKeyF3(R,key) : key in keys ];
        cusp_count := #triples;
        triples cat:= [ [R | item[1],0,0] : item in EF3_QuarticOrbitData(R) | item[1] ne 0 ];
        models := [ [R | 0,W[1],0,W[2],W[3]] : W in triples ];
    end if;
    return models,cusp_count;
end function;

function reference_singular_cells(models)
    cells := [0,0,0,0];
    assert not IsEmpty(models);
    for W in models do
        kind,minimum := reference_properties(W);
        assert kind ne "smooth";
        i := kind eq "regular cusp" select 1 else 3;
        if minimum eq "nonminimal somewhere" then i +:= 1; end if;
        cells[i] +:= 1;
    end for;
    return cells;
end function;
