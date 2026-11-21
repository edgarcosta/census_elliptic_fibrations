function torsion_samples(F)
    t := F.1;
    p := Characteristic(F);
    samples := p eq 2 select [
        [F | 0, 0, 1, t^6, 0],
        [F | 0, 0, 1, t^7, 0],
        [F | t, 0, t^6+1, 0, 0],
        [F | t^2, 0, 1, t^8, t^6],
        [F | 1, 0, t^6, 1, t^6],
        [F | t^2, 0, 1, 0, 0]] else [
        [F | 0, 0, 0, t^8, t+1],
        [F | 0, t^3+t, 0, t^8+1, 0],
        [F | 0, t^4, 0, 0, 2],
        [F | 0, t^4+1, 0, t^4, 0],
        [F | 1, t^4, t, t^4+t, -t^2],
        [F | 0, 1+t^4, 0, t^2-t^4, t^4]];
    R<s, u> := PolynomialRing(GF(p), 2);
    keys := p eq 2 select [1346994688] else [
        879831930519, 65172735594, 2542281946632, 75547326549];
    for key in keys do
        Append(~samples, aInvariants(EllipticCurve(aInvariantsFromKey(R, key), F)));
    end for;
    return samples;
end function;
