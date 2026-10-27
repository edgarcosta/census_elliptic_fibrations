SetQuitOnError(true);
AttachSpec("spec");
import "src/census.m": prefixes, fibre, equation_weight;
load "reference/support/long_models.m";
// Public row totals pin each raw-prefix multiplicity. Engine dimensions and
// stored stabilizers also agree with independently enumerated translation images.
procedure test_public_equation_weights()
    for p in [2,3] do
     R<s,t> := PolynomialRing(GF(p),2);
     PP := FibrationPrefixes(R);
     private := prefixes(R);
     assert [ <P[1],P[2]> : P in private ] eq PP;
     _,_,rows := NumberOfEquations(R);
     base := [ Matrix(g) : g in GL(2,GF(p)) ];
     shifts := math_forms(R,p eq 2 select 2 else 4);
     sextics := p eq 2 select math_forms(R,6) else [R|];
     raw := 0;
     for i in [1..#PP] do
      a,b := Explode(PP[i]);
      stabilizer := 0;
      for g in base do
       coords := [g[1,1]*s+g[1,2]*t,g[2,1]*s+g[2,2]*t];
       ag := Evaluate(a,coords);
       if ag ne a then continue; end if;
       bg := Evaluate(b,coords);
       if p eq 2 then
        stabilizer +:= #[ r : r in shifts | bg+a*(r^2+a*r) eq b ];
       elif a eq 0 then
        if bg eq b then stabilizer +:= 1; end if;
       elif IsDivisibleBy(bg-b,a) then
        stabilizer +:= 1;
       end if;
      end for;
      assert stabilizer eq #private[i][3];
      if p eq 2 then
       image := { <a*v,v^2+b*v> : v in sextics };
       mult := 48 div stabilizer;
      else
       image := a eq 0 select { <r^3+b*r> : r in shifts } else { <R!0> };
       mult := (48 div stabilizer)*(a eq 0 select 1 else 243);
      end if;
      F := fibre(R,private[i]);
      assert #image eq p^((p eq 2 select 22 else 13)-F`Dimension);
      assert equation_weight(R,private[i],F) eq mult*#image;
      assert &+[ rows[i][2][cell] : cell in Keys(rows[i][2]) ] eq
          mult*p^(p eq 2 select 22 else 13);
      raw +:= mult;
     end for;
     assert raw eq (p eq 2 select 2^10 else 3^14);
    end for;
end procedure;

try
    test_public_equation_weights();
catch e
    printf "MISMATCH: public_equation_weights %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
