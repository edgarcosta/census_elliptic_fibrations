function fixed_models(R)
    p := Characteristic(R);
    F<t> := FunctionField(GF(p));
    s := R.1;
    u := R.2;
    weights := [2, 4, 6, 8, 12];
    equations := p eq 2 select [
        [F | t, 0, t^6+1, 0, 0],
        [F | t, 0, 1, 0, t^12],
        [F | 0, 0, 1, t^6, 0],
        [F | 0, 0, 1, t^7, 0]] else [
        [F | 0, t^3+t, 0, t^8+1, 0],
        [F | 0, t^4, 0, 1, t^12],
        [F | 0, 0, 0, t^8, t+1]];
    models := [[R | &+[R | Coefficient(Numerator(a[i]), j)*s^j*u^(weights[i]-j)
        : j in [0..Degree(Numerator(a[i]))]] : i in [1..5]] : a in equations];
    keys := p eq 2 select [4884480, 11821184] else [372699133912, 320009441124];
    models cat:= [aInvariantsFromKey(R, key) : key in keys];
    return models;
end function;
