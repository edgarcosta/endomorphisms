// Regression cases for the upperbounds Magma port: one case per bug we have
// actually made and fixed, so that bug cannot come back. Every case here must
// fail on the pre-fix code; one that passes both ways belongs in another suite.
// New bugs get new cases here, each naming the bug it guards.

// Requires PARI/gp on PATH: the center-field detection (FieldIntersection ->
// SubfieldsPolynomials -> Polredabs) needs it. If gp is missing we print a
// notice and exit rather than report spurious failures.

// Silent on success. Run via:
//   ./tests/run.sh regressions (from the repository root)

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

if #Pipe("command -v gp || true", "") eq 0 then
    print "Test-UpperBounds-Regressions: gp (PARI/gp) not on PATH; skipping.";
    exit;
end if;

// ----- Bug endomorphisms-6qm (fixed in 88e7375), 1 of 2 -----
// The escalation returned "unsound" and stopped climbing at the first rung whose
// bound missed the target. The factor and centre bounds are conditional on eta
// and t (UpperBounds.m:167), so a low rung can legitimately miss the truth.

// 10368.ba1, truth [CC, CC]: B = 10 has only 2 surviving primes and yields
// [M_2(RR), RR], which misses [CC, CC]; B = 20 has 6 primes and is sharp.
// Pre-fix this returned "unsound" at B = 10 and never reached B = 20.
C10368 := HyperellipticCurve(R ! [0,0,3,0,9,0,4], R ! [1]);
bound, B, status := RealRepresentationBound(C10368, ["CC", "CC"]);
error if status ne "sharp" or B ne 20,
    Sprintf("10368.ba1: expected sharp at B = 20, got %o at B = %o with bound %o",
            status, B, bound);

// ----- Bug endomorphisms-6qm (fixed in 88e7375), 2 of 2 -----
// Counterpart of the fix above: with an unsound rung made non-fatal, "unsound"
// must still fire at the top of the ladder, where the ~130 usable primes below
// 800 make a non-containment a real defect. Guards it against never firing.

// 256.a.512.1, truly M_2(RR), fed the false target [CC, RR]: C x R has no unital
// injection into M_2(R), as the two central idempotents would be complementary
// and rank one, forcing the C summand into P M_2(R) P = R. No rung can match, so
// the climb must run the whole ladder and report unsound at its top rung, 800.
C256 := HyperellipticCurve(R ! [0,-1,1,1,-3,2], R ! [1]);
bound, B, status := RealRepresentationBound(C256, ["CC", "RR"]);
error if status ne "unsound" or B ne 800,
    Sprintf("256.a.512.1 with false target [CC, RR]: expected unsound at B = 800, got %o at B = %o with bound %o",
            status, B, bound);

// The same curve with its true target still resolves, so the case above pins a
// false target being caught, not the curve being unreachable.
bound, B, status := RealRepresentationBound(C256, ["M_2(RR)"]);
error if status ne "sharp",
    Sprintf("256.a.512.1 with true target [M_2(RR)]: expected sharp, got %o at B = %o with bound %o",
            status, B, bound);

// ----- Bug endomorphisms-91i (fixed in 12e59e6) -----
// A result was called a valid upper bound whenever its total real dimension
// exceeded the truth's, but containment needs a unital injective embedding.
// RealRepresentationEmbeds in upperbounds/match.m is the Magma side of it.

// Rows 1 and 4 are what a dimension-only test gets wrong: [CC, RR] has dimension
// 3 < 4 yet does not embed in M_2(RR), and [HH] matches the dimension exactly
// without being isomorphic to it.
embeds := [*
    <["CC", "RR"],       ["M_2(RR)"], false>,
    <["RR", "RR"],       ["M_2(RR)"], true>,
    <["RR", "RR", "RR"], ["M_2(RR)"], false>,
    <["HH"],             ["M_2(RR)"], false>,
    <["CC", "CC"],       ["M_2(CC)"], true>,
    <["M_2(RR)"],        ["M_2(CC)"], true>
*];

for entry in embeds do
    target, bnd, expected := Explode(entry);
    got := RealRepresentationEmbeds(target, bnd);
    error if got ne expected,
        Sprintf("RealRepresentationEmbeds(%o, %o): expected %o, got %o",
                target, bnd, expected, got);
end for;

// ----- Ladder capped at B = 200 (fixed by appending the 400 and 800 rungs) -----
// Only three of 3,420,837 curves exhausted the old ladder: each returns [CC, RR]
// at B = 200, sound since R x R embeds in C x R, but not sharp. One centre
// appears imaginary quadratic instead of Q; all three become sharp at 400.

exhausted200 := [*
    <"177674.e1",  [0,49,-37,143,47,81],        [0,1,1]>,
    <"203056.bc1", [-381,388,-58,464,139,44,8], []>,
    <"406112.o1",  [95,-97,14,-116,-35,-11,-2], [1,0,1]>
*];

for entry in exhausted200 do
    label, fc, hc := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound, B, status := RealRepresentationBound(C, ["RR", "RR"]);
    error if status ne "sharp" or B ne 400,
        Sprintf("%o: expected sharp at B = 400, got %o at B = %o with bound %o",
                label, status, B, bound);
end for;

// ----- Coverage of the B = 200 tail, NOT regressions -----
// Measured hard cases; only 0.09 percent of the corpus reaches B = 200. Each
// asserts sharp and a rung no larger than the one recorded, so a change needing
// more primes fails here while an improvement needing fewer still passes.

// <label, f-coeffs, h-coeffs, truth tokens, max expected rung>.
tail := [*
    <"1000377.e1",  [0,16,88,-39,-264,144], [],        ["RR","RR"], 200>,
    <"136500.bc4",  [312523400,910850700,-913084905,-123883730,-206492181,-180183612,3299858], [0,0,1,1], ["RR","RR"], 200>,
    <"245700.ee6",  [-32426971668000,0,-1041100665919,0,31878491585,0,89504640], [0,1], ["RR","RR"], 200>,
    <"362624.c1",   [-28,60,11,-80,21,20,-4], [],      ["RR","RR"], 200>,
    <"142002.m1",   [-96312,-69027,119,-1676,17,-9], [0,1,1], ["CC","RR"], 200>,
    <"199969.b1",   [378,378,693,392,339,92,43], [],   ["CC","RR"], 200>,
    <"486717.d1",   [-183425,73312,-86101,25984,-14007,1218,-406], [], ["CC","RR"], 200>,
    // Intersection-starvation witness: at B = 50 only p = 31 and 37 reach
    // eta_lower and both carry the discriminant pair {-11,-3}, so the
    // intersection is a no-op; p = 59 at B = 100 breaks it and gives [CC, RR].
    <"1000912.l1",  [82,0,203,0,77,0,-77], [1],        ["CC","RR"], 100>,
    <"259840.bm3",  [-255,1774,-1271,156,-106,-28], [], ["CC","RR"], 100>,
    <"153664.cq5",  [168,0,126,0,693,0,168], [],       ["CC","CC"], 50>,
    <"165888.f1",   [-1,-5,6,44,6,-5,-1], [1,1,1,1],   ["CC","CC"], 20>,
    <"20736.bm3",   [2,0,-6,0,12,0,-8], [],            ["CC","CC"], 20>
*];

for entry in tail do
    label, fc, hc, truth, maxB := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound, B, status := RealRepresentationBound(C, truth);
    error if status ne "sharp",
        Sprintf("%o: expected sharp, got %o at B = %o with bound %o, truth %o",
                label, status, B, bound, Sort(truth));
    error if B gt maxB,
        Sprintf("%o: sharp at B = %o, above the recorded max %o (bound %o)",
                label, B, maxB, bound);
end for;

// ----- Coverage of the field intersection, NOT regressions -----
// Row count does not predict cost: 850950.n1 has the corpus maximum 24 rows and
// is the cheapest here, while 840693.b1 has 11 and is 5x slower. Cost tracks the
// rows the walk visits before the running intersection collapses to Q.
coverageVerified := [*
    <"840693.b1",  [-40334,-100506,434601,-314128,-15833,-220,-1], [],        200, ["CC","RR"]>,
    <"486717.d1",  [-183425,73312,-86101,25984,-14007,1218,-406], [],        200, ["CC","RR"]>,
    <"850950.n1",  [-158720,0,-37871,0,-2717,0,-61], [0,0,1],               200, ["RR","RR"]>,
    <"121680.el3", [0,-39,97,-30,-28,-1], [0,1,0,1],                         200, ["RR","RR"]>,
    <"409600.bh1", [0,2,0,1,-1,-1,-1], [],                                  50, ["RR","RR"]>,
    <"422500.v1",  [20,-90,-434,-218,316,210,20], [0,1,1],                   50, ["M_2(RR)"]>
*];

for entry in coverageVerified do
    label, fc, hc, maxB, expected := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp",
        Sprintf("%o: expected sharp, got %o at B = %o with bound %o, truth %o",
                label, status, B, bound, Sort(expected));
    error if B gt maxB,
        Sprintf("%o: sharp at B = %o, above the recorded max %o (bound %o)",
                label, B, maxB, bound);
end for;
