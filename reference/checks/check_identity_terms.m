SetQuitOnError(true);
AttachSpec("spec");
import "src/census.m": prefixes, fibre;
import "src/engine.m": fibre_counts, fibre_equation_counts;

// Boolean locus patterns are an engine contract, beyond the public six cells.
// check_equation_rows.m also checks these fibres through NumberOfEquations.
procedure test_identity_locus_patterns()
    for p in [2,3] do
     R := PolynomialRing(GF(p),2);
     PP := prefixes(R);
     indices := p eq 2 select [1,2,33] else [1,3,1330];
     for i in indices do
      F := fibre(R,PP[i]);
      order := #F`Group;
      api := fibre_equation_counts(F);
      assert #F`Group eq order;
      F`Group := sub< F`Group | >;
      direct := AssociativeArray();
      for flags in CartesianPower([false,true],3) do direct[<b : b in flags>] := 0; end for;
      for v in VectorSpace(GF(p),F`Dimension) do
       pattern := < &or[ (v-A[1]) in A[2] : A in locus ] : locus in F`Loci >;
       direct[pattern] +:= 1;
      end for;
      actual := fibre_counts(F);
      assert Keys(actual) eq Keys(direct) and Keys(api) eq Keys(direct);
      for pattern in Keys(direct) do
       assert actual[pattern] eq direct[pattern];
       assert api[pattern] eq direct[pattern];
      end for;
     end for;
    end for;
end procedure;

try
    test_identity_locus_patterns();
catch e
    printf "MISMATCH: identity_locus_patterns %o\n",e;
    quit 1;
end try;
if assigned resultpath then Write(resultpath,"REFERENCE_CHECK_OK"); end if;
quit 0;
