import "tests/support/_archive.m": test_archive;
import "tests/support/_archive_oracle.m": test_archive_canonicalization;

// Two complete F_3 prefixes exercise both canonicalization branches.
procedure test_archive_invariants_f3()
    test_archive_canonicalization(3);
    expected := [
        <"all", 1594638, "37f8912550459978713f03631cf592873fbab915493b388d6811c2dc79479f0e">,
        <"smooth_minimal", 1594323, "f9ce481ae63ecfb0103ee69646d9d6a51bb8451cbd61010a78d4f6ac99b5094a">,
        <"smooth_nonminimal", 0, "f879835f97d705295a6ec2043ceccacdab5f27f792fd342313d01e71eba7c2d2">,
        <"delta_zero", 315, "c250c19f0e4e7f9fead833b03255ddde6999b73183b0eb181406b6f00ab990e2">
    ];
    test_archive(3, {0, 21870}, expected);
end procedure;

test_archive_invariants_f3();
