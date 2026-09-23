// Coverage from the finished 3,420,837-curve LMFDB validation run.
// Curves partition into (stratum, B-rung) cells: stratum is geom_end_alg,
// and B is the rung where the escalating bound matched the truth. Of these,
// 27 cells are populated; each contributes its smallest-conductor curve.

// The listed curves are the existing cell witnesses, rebased to their exact
// first sharp rung on the current schedule.

// Requires PARI/gp on PATH for center-field detection via Polredabs.
// Silent on success. Run from the root via: ./tests/run.sh corpus-cells

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

if #Pipe("command -v gp || true", "") eq 0 then
    print "corpus-cells: gp (PARI/gp) not on PATH; skipping.";
    exit;
end if;

// <label, stratum, recorded B, expected tokens, f-coeffs, h-coeffs (low to high)>.
cells := [*
    <"28561.c1", "CM", 8, ["CC", "CC"], [0, 0, -2, 2, 2, -3, -2], [1]>,
    <"3125.a1", "CM", 16, ["CC", "CC"], [0, 0, 0, 0, 0, -1], [1]>,
    <"20736.ce1", "QM", 8, ["M_2(RR)"], [0, -2, 0, 5, -3, -5, 1], [1, 1, 1, 1]>,
    <"8100.o1", "QM", 16, ["M_2(RR)"], [1, 3, -42, 43, 21, -60, -28], [1]>,
    <"8192.d1", "CM x CM", 32, ["CC", "CC"], [-2, 0, 8, 0, -9, 0, 2], []>,
    <"76832.bb1", "CM x CM", 32, ["CC", "CC"], [0, -7, -30, 15, 9, -2, -1], [0, 1, 0, 1]>,
    <"729.a1", "M_2(CM)", 8, ["M_2(CC)"], [16, 0, 0, -5], [0, 0, 0, 1]>,
    <"576.a1", "M_2(CM)", 16, ["M_2(CC)"], [-2, 0, 0, 5], [0, 0, 0, 1]>,
    <"290521.b1", "M_2(CM)", 32, ["M_2(CC)"], [3993, 0, 1361, 0, -495, 0, -192], [0, 1]>,
    <"529.a1", "RM", 8, ["RR", "RR"], [0, 1, 0, -1, 0, -1], [1, 0, 1, 1]>,
    <"1521.a1", "RM", 32, ["RR", "RR"], [2, -3, 4, -2, 1], [1, 0, 0, 1]>,
    <"8281.b1", "RM", 32, ["RR", "RR"], [-2, -4, -5, -4, -4, -1, -1], [1, 1, 1, 1]>,
    <"7225.e1", "RM", 64, ["RR", "RR"], [14, -36, -32, 65, 16, 23, 14], [1, 1, 1]>,
    <"104544.hh1", "RM", 128, ["RR", "RR"], [9, -18, -24, 45, 12, -28, 5], [0, 0, 1]>,
    <"121.a1", "M_2(Q)", 8, ["M_2(RR)"], [-2, 4, 2, 5, 2, 1], [1, 1, 1, 1]>,
    <"196.a1", "M_2(Q)", 16, ["M_2(RR)"], [1, 3, 6, 7, 6, 3, 1], [0, 1, 1]>,
    <"16900.f1", "M_2(Q)", 32, ["M_2(RR)"], [4, 42, 63, -44, -87, -18, 4], [0, 1, 1]>,
    <"6859.a1", "CM x Q", 8, ["CC", "RR"], [0, -1, -1, -11, -4, -16, 16], [1]>,
    <"378.a1", "CM x Q", 32, ["CC", "RR"], [-13896, -1861, -1185, -38, -23, 2], [1, 0, 1]>,
    <"448.a1", "CM x Q", 32, ["CC", "RR"], [-7, 0, 0, 0, 1], [0, 1, 0, 1]>,
    <"1026.d1", "CM x Q", 128, ["CC", "RR"], [-307, -675, -789, -1240, -338, -490, 118], [0, 1, 0, 1]>,
    <"47334.m1", "CM x Q", 128, ["CC", "RR"], [-24005, 7583, 7735, 101, -457, -609, -203], [1, 1, 1]>,
    <"363.a1", "Q x Q", 8, ["RR", "RR"], [0, -1, 0, -2, 4, -2], [0, 1, 0, 1]>,
    <"294.a1", "Q x Q", 32, ["RR", "RR"], [0, 0, 1, -1, 1], [1, 0, 0, 1]>,
    <"255.a1", "Q x Q", 64, ["RR", "RR"], [0, 16, 1, 8, 0, 1], [0, 0, 1]>,
    <"510.a1", "Q x Q", 128, ["RR", "RR"], [-58, 108, -62, -214, 404, -357, 119], [1, 1, 1]>,
    <"2080.c1", "Q x Q", 128, ["RR", "RR"], [-13, 4, -25, 57, -117, 36, -3], [1, 1]>
*];

for entry in cells do
    label, stratum, maxB, expected, fc, hc := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B ne maxB,
        Sprintf("%o (%o): expected %o, sharp at B = %o; got %o at B = %o with bound %o",
                label, stratum, Sort(expected), maxB, status, B, bound);
end for;
