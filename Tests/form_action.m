import "src/forms.m": form_action, form_from_vector, base_change;

// Pins the SymmetricPower row action against direct base substitution.
procedure test_form_action()
    for p in [2,3] do
        k := GF(p);
        R<s,t> := PolynomialRing(k,2);
        G := GL(2,k);
        for d in [2,4,6,8,12] do
            rho := form_action(R,G,d);
            V := VectorSpace(k,d+1);
            for g in G do
                for i in [1..d+1] do
                    f := form_from_vector(R,d,V.i*rho(g));
                    assert f eq base_change(s^(d-i+1)*t^(i-1),Matrix(g));
                end for;
            end for;
        end for;
    end for;
end procedure;

test_form_action();
