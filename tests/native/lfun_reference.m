SetQuitOnError(true);
SetNthreads(1);
p := StringToInteger(p);
F<t> := FunctionField(GF(p));
a := [F | PolynomialRing(GF(p))!IntegerToSequence(StringToInteger(s), p)
    : s in Split(Read(inputpath), " \n")];
L := LFunction(EllipticCurve(a));
Write(outputpath, Join([Sprint(Integers()!c) : c in Coefficients(L)], " "));
quit;
