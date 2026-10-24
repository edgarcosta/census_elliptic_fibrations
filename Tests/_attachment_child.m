SetQuitOnError(true);
AttachSpec(specpath);
for p in [2,3] do
    R := PolynomialRing(GF(p),2);
    assert #FibrationPrefixes(R) eq (p eq 2 select 88 else 1860);
end for;
Write(resultpath,"ATTACHMENT_OK");
exit 0;
