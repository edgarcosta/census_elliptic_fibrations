function torsion_samples(F)
    t := F.1;
    p := Characteristic(F);
    samples := p eq 2 select [
        <[F | 0, 0, 1, t^6, 0], []>,
        <[F | 0, 0, 1, t^7, 0], []>,
        <[F | t, 0, t^6+1, 0, 0], [3]>,
        <[F | t^2, 0, 1, t^8, t^6], [2]>,
        <[F | 1, 0, t^6, 1, t^6], [2]>,
        <[F | t^2, 0, 1, 0, 0], [6]>] else [
        <[F | 0, 0, 0, t^8, t+1], []>,
        <[F | 0, t^3+t, 0, t^8+1, 0], [2]>,
        <[F | 0, t^4, 0, 0, 2], [3]>,
        <[F | 0, t^4+1, 0, t^4, 0], [2, 4]>,
        <[F | 1, t^4, t, t^4+t, -t^2], [2, 4]>];
    R<s, u> := PolynomialRing(GF(p), 2);
    key := p eq 2 select 1346994688 else 75547326549;
    Append(~samples, <[F | Evaluate(a, [t, F!1])
        : a in aInvariantsFromKey(R, key)], []>);
    return samples;
end function;
