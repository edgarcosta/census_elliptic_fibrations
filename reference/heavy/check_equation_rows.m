SetQuitOnError(true);
AttachSpec("spec");
load "reference/support/long_models.m";
load "reference/support/heavy.m";
load "reference/support/parallel.m";

// Enumerates quotient equations with independently constructed translation images.
procedure test_equation_row(p,index)
    R<s,t> := PolynomialRing(GF(p),2);
    k := BaseRing(R);
    P := FibrationPrefixes(R)[index];
    first,second := Explode(P);
    raw_prefixes := {};
    shifts := math_forms(R,p eq 2 select 2 else 4);
    for g in GL(2,k) do
     a := Evaluate(first,[g[1,1]*s+g[1,2]*t,g[2,1]*s+g[2,2]*t]);
     b := Evaluate(second,[g[1,1]*s+g[1,2]*t,g[2,1]*s+g[2,2]*t]);
     for r in shifts do
      if p eq 2 then Include(~raw_prefixes,<a,b+a*(r^2+a*r)>);
      else Include(~raw_prefixes,<a,b-a*r>); end if;
     end for;
    end for;
    V := VectorSpace(k,p eq 2 select 22 else 13);
    if p eq 2 then
     images := { V!(math_coefficients(first*v,8) cat math_coefficients(v^2+second*v,12)) : v in math_forms(R,6) };
    elif first eq 0 then
     images := { V!math_coefficients(r^3+second*r,12) : r in math_forms(R,4) };
    else
     images := { V!0 };
    end if;
    S := sub< V | Setseq(images) >;
    assert #S eq #images;
    U := sub< V | [ ReduceVector(S,v) : v in Basis(V) ] >;
    assert Dimension(S meet U) eq 0;
    assert Dimension(S+U) eq Dimension(V);
    assert #U*#images eq #V;
    cells := [ <a,b> : a in ["smooth","regular cusp","nonregular"],
     b in ["everywhere minimal","nonminimal somewhere"] ];
    counts := [0 : c in cells];
    checked := 0;
    for v in U do
     z := p eq 2 select [ math_form(R,Eltseq(v)[1..9],8),
         math_form(R,Eltseq(v)[10..22],12) ] else [ math_form(R,Eltseq(v),12) ];
     W := p eq 2 select [R | first,0,second,z[1],z[2]] else [R | 0,first,0,second,z[1]];
     a,b := FibrationProperties(W);
     counts[Index(cells,<a,b>)] +:= #images;
     checked +:= 1;
    end for;
    assert &+counts eq #V;
    counts := [ n*#raw_prefixes : n in counts ];
    total,table,rows := NumberOfEquations(R : PrefixIndices := [index]);
    assert #rows eq 1 and rows[1][1] eq index;
    assert total eq &+counts;
    for j in [1..#cells] do
     assert table[cells[j]] eq counts[j];
     assert rows[1][2][cells[j]] eq counts[j];
    end for;
end procedure;

try
    processes := Min(heavy_processes(),6);
    value := GetEnv("REFERENCE_SHARD");
    if value eq "" and processes gt 1 then
        reference_parallel_check("reference/heavy/check_equation_rows.m",processes);
    else
        shard := value eq "" select 1 else heavy_shard(processes);
        targets := [<2,1>,<2,2>,<2,33>,<3,1>,<3,3>,<3,1330>];
        for i in [1..#targets] do
            if (i-1) mod processes eq shard-1 then
                test_equation_row(targets[i][1],targets[i][2]);
            end if;
        end for;
    end if;
catch e
    printf "MISMATCH: equation row enumeration %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
