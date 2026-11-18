SetQuitOnError(true);
SetNthreads(1);
AttachSpec("spec");
import "tests/support/_native_derived_fields.m": derived_fields;
import "tests/support/magma/_derived_fixtures.m": gf4096_fixture, gf4_fixture;

procedure write_constant_field_controls(path)
    file := Open(path, "w");
    delete file;
    E, L, rows := gf4096_fixture();
    result := FibrationInvariants(E, L, rows, []);
    assert LeadingCoefficient(result`P2) eq 2^264;
    assert result`PointCounts[5] eq 1329227995784915898268080161630978049;
    fprintf path, "%o\n", Join(derived_fields(result), ",");
    E, L, rows := gf4_fixture();
    result := FibrationInvariants(E, L, rows, []);
    assert result`BrauerHeight eq 3;
    fprintf path, "%o\n", Join(derived_fields(result), ",");
end procedure;

write_constant_field_controls(referencepath);
quit;
