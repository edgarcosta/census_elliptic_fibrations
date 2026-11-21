import "tests/support/_archive.m": test_archive;
import "tests/support/_archive_oracle.m": test_archive_canonicalization;

// All F_2 prefixes compare complete orbit multisets and archived classification.
procedure test_archive_fixtures_f2()
    test_archive_canonicalization(2);
    expected := [
        <"all", 1402468, "58faf27cefa9467bde09429c7dd0c8f723b06cd77c3cd582c2f522859570d5fa">,
        <"smooth_minimal", 1397598, "f0449084059f1fccf354537604cbb416f98a62efdc05a99a3b21394e66b5d964">,
        <"smooth_nonminimal", 4134, "c4af195e25a1a05e391648e47540d4f3e96b639d6a5ad90f52d1614306539b62">,
        <"delta_zero", 736, "886de781dc63795a76c0df1ee857f48bd4f91c33eb28f34238d751d83697d175">
    ];
    test_archive(2, {Integers() | }, expected);
end procedure;

test_archive_fixtures_f2();
