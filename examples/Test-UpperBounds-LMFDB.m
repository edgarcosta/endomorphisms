// LMFDB genus-2 corpus for the upperbounds Magma port.
//
// Curves pulled from the LMFDB (g2c_curves joined to g2c_endomorphisms), a
// small-discriminant batch per geometric real-endomorphism-algebra stratum.
// The LMFDB minimal model eqn = [[f-coeffs], [h-coeffs]] is low->high, for
// y^2 + h(x) y = f(x), and decodes as HyperellipticCurve(f, h).
//
// Each row hands its factorsRR_geom to the escalating overload
//   RealRepresentationBound(C::Crv, target::SeqEnum) -> bound, B, status
// which climbs B over 10, 20, 50, 100, 200 and stops at the first match, so
// there is no per-row B any more; we assert it reports "sharp".
//
// That pins two things at once: the bound is reached, and it is never
// contradicted. Status "unsound" means the truth does not embed in the
// computed bound, which the theorem forbids, so it stops the climb.
//
// Requires PARI/gp on PATH: the center-field detection (FieldIntersection ->
// SubfieldsPolynomials -> Polredabs) needs it, and without it the CM strata
// collapse (e.g. ["RR","CC"] -> ["RR","RR"]). If gp is not found we print a
// notice and skip rather than report spurious failures.
//
// Sage ships a usable gp; e.g. run with
//   PATH="$(dirname $(sage -sh -c 'command -v gp')):$PATH".
//
// Silent on success. Set ENDO_TEST_VERBOSE to a non-empty value to print the
// minimal sharp B per curve.
//
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds-LMFDB.m Utils/Exit.m
// (kept separate from Test-UpperBounds.m, whose helper assertions pin
//  non-polredabs polynomial forms and so are written for the gp-free path.)

AttachSpec("../endomorphisms/magma/spec");
AttachSpec("/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec");

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

if #Pipe("command -v gp || true", "") eq 0 then
    print "Test-UpperBounds-LMFDB: gp (PARI/gp) not on PATH; skipping corpus.";
    exit;
end if;

verbose := GetEnv("ENDO_TEST_VERBOSE") ne "";

// <label, f-coeffs, h-coeffs, factorsRR_geom (flat)>.
corpus := [*
    <"249.a.249.1",     [0,1,1],          [1,0,0,1], ["RR"]>,
    <"277.a.277.1",     [0,-1,-1],        [1,1,1,1], ["RR"]>,
    <"295.a.295.1",     [0,0,-1],         [1,0,0,1], ["RR"]>,
    <"529.a.529.1",     [0,0,0,0,0,-1],   [1,1,0,1], ["RR","RR"]>,  // RM, geom-simple
    <"841.a.841.1",     [2,1,3,1,1],      [0,1,1,1], ["RR","RR"]>,  // RM, geom-simple
    <"294.a.294.1",     [0,0,1,0,1],      [1,0,0,1], ["RR","RR"]>,  // E1 x E2 split
    <"847.a.847.1",     [0,0,1,1,1],      [1,1,1,1], ["RR","RR"]>,  // E1 x E2 split
    <"448.a.448.1",     [-7,0,0,0,1],     [0,1,0,1], ["RR","CC"]>,  // E x E_CM
    <"1331.a.1331.1",   [1,3,2,-1,-1],    [0,0,0,1], ["RR","CC"]>,  // E x E_CM
    <"256.a.512.1",     [0,-1,1,1,-3,2],  [1],       ["M_2(RR)"]>,  // E^2, non-CM
    <"3721.a.3721.1",   [0,1,3,1,-1],     [1,1,0,1], ["M_2(RR)"]>,  // E^2, non-CM
    <"11664.a.11664.1", [0,0,0,0,0,0,-1], [1],       ["M_2(CC)"]>,  // E^2, CM
    <"4096.b.65536.1",  [0,-1,0,0,0,1],   [],        ["M_2(CC)"]>   // E^2, CM
*];

for entry in corpus do
    label, fc, hc, expected := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound, B, status := RealRepresentationBound(C, expected);
    error if status ne "sharp",
        Sprintf("%o: status %o at B = %o, expected %o, got %o",
                label, status, B, Sort(expected), bound);
    if verbose then
        printf "%-18o sharp at B = %3o  %o\n", label, B, bound;
    end if;
end for;

// Geometrically simple quartic CM: the centre bound must name the CM field, not
// its real quadratic subfield, which would report [RR, RR] and halve the bound.
// Kept on the plain overload, as a pin on the fix rather than on escalation.
C3125 := HyperellipticCurve(x^5, R ! 1);  // 3125.a.3125.1: y^2 + y = x^5, CM by Q(zeta_5)
error if Sort(&cat RealRepresentationBound(C3125, 200)) ne ["CC", "CC"],
    "3125.a.3125.1: LMFDB factorsRR_geom is [CC, CC]";
