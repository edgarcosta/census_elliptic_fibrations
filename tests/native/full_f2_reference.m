SetQuitOnError(true);
SetNthreads(1);
SetColumns(0);
AttachSpec("spec");
import "tests/_native_fixed_models.m": fixed_models;

file := Open(keypath, "w");
Puts(file, "key,cell,prefix");
R<s, t> := PolynomialRing(GF(2), 2);
kinds := ["smooth", "regular cusp", "nonregular"];
minimalities := ["everywhere minimal", "nonminimal somewhere"];
counts := [0 : i in [1..6]];
for i in [1..#FibrationPrefixes(R)] do
    for k in [1..3] do
        for m in [1..2] do
            cell := 2*(k-1) + m-1;
            keys := FibrationRepresentativeKeys(R, i :
                GenericFibreKinds := [kinds[k]], Minimalities := [minimalities[m]]);
            for key in keys do
                Puts(file, Sprintf("%o,%o,%o", key, cell, i));
            end for;
            counts[cell+1] +:= #keys;
        end for;
    end for;
end for;
Flush(file);
delete file;
assert counts eq [1397598, 4134, 680, 38, 7, 11];

file := Open(referencepath, "w");
Puts(file, FibrationInvariantCSV());
R<s, t> := PolynomialRing(GF(3), 2);
for W in fixed_models(R) do
    Puts(file, FibrationInvariantCSV(FibrationInvariantRow(FibrationKey(W), W)));
end for;
Flush(file);
quit;
