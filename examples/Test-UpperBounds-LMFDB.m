// LMFDB genus-2 corpus for the upperbounds Magma port.
//
// Curves pulled from the LMFDB (g2c_curves joined to g2c_endomorphisms), a
// small-discriminant batch per geometric real-endomorphism-algebra stratum.
// For each curve we decode the LMFDB minimal model
//   eqn = [[f-coeffs], [h-coeffs]]   (low->high; the model is y^2 + h(x) y = f(x))
// into HyperellipticCurve(f, h) and assert that the flattened upper bound matches
// LMFDB's factorsRR_geom exactly. RealRepresentationBound returns one list per
// simple geometric factor, so its concatenation is the flat multiset that LMFDB
// stores; we compare as sorted multisets:
//   Sort(&cat RealRepresentationBound(C, B)) eq Sort(factorsRR_geom).
//
// Requires PARI/gp on PATH: the center-field detection (FieldIntersection ->
// SubfieldsPolynomials -> Polredabs) needs it, and without it the CM strata
// collapse (e.g. ["RR","CC"] -> ["RR","RR"]). If gp is not found we print a
// notice and skip rather than report spurious failures. Sage ships a usable gp;
// e.g. run with  PATH="$(dirname $(sage -sh -c 'command -v gp')):$PATH".
//
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds-LMFDB.m Utils/Exit.m
// (kept separate from Test-UpperBounds.m, whose helper assertions pin
//  non-polredabs polynomial forms and so are written for the gp-free path.)

AttachSpec("../endomorphisms/magma/spec");
AttachSpec("/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec");

R<x> := PolynomialRing(Rationals());

if #Pipe("command -v gp || true", "") eq 0 then
    print "Test-UpperBounds-LMFDB: gp (PARI/gp) not on PATH; skipping corpus.";
    exit;
end if;

// <label, f-coeffs, h-coeffs, factorsRR_geom (flat), B>.  Sharp at B = 200.
corpus := [*
    <"249.a.249.1",     [0,1,1],          [1,0,0,1], ["RR"],      200>,
    <"277.a.277.1",     [0,-1,-1],        [1,1,1,1], ["RR"],      200>,
    <"295.a.295.1",     [0,0,-1],         [1,0,0,1], ["RR"],      200>,
    <"529.a.529.1",     [0,0,0,0,0,-1],   [1,1,0,1], ["RR","RR"], 200>,  // RM, geom-simple
    <"841.a.841.1",     [2,1,3,1,1],      [0,1,1,1], ["RR","RR"], 200>,  // RM, geom-simple
    <"294.a.294.1",     [0,0,1,0,1],      [1,0,0,1], ["RR","RR"], 200>,  // E1 x E2 split
    <"847.a.847.1",     [0,0,1,1,1],      [1,1,1,1], ["RR","RR"], 200>,  // E1 x E2 split
    <"448.a.448.1",     [-7,0,0,0,1],     [0,1,0,1], ["RR","CC"], 200>,  // E x E_CM
    <"1331.a.1331.1",   [1,3,2,-1,-1],    [0,0,0,1], ["RR","CC"], 200>,  // E x E_CM
    <"256.a.512.1",     [0,-1,1,1,-3,2],  [1],       ["M_2(RR)"], 200>,  // E^2, non-CM
    <"3721.a.3721.1",   [0,1,3,1,-1],     [1,1,0,1], ["M_2(RR)"], 200>,  // E^2, non-CM
    <"11664.a.11664.1", [0,0,0,0,0,0,-1], [1],       ["M_2(CC)"], 200>,  // E^2, CM
    <"4096.b.65536.1",  [0,-1,0,0,0,1],   [],        ["M_2(CC)"], 200>   // E^2, CM
*];

for entry in corpus do
    label, fc, hc, expected, B := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    bound := RealRepresentationBound(C, B);
    // Non-empty => the putative upper bound was established for this curve.
    error if #bound eq 0, Sprintf("%o: no upper bound established", label);
    error if Sort(&cat bound) ne Sort(expected),
        Sprintf("%o: expected %o, got %o", label, Sort(expected), Sort(&cat bound));
end for;

// Geometrically simple quartic CM: the centre bound must name the CM field, not
// its real quadratic subfield, which would report [RR, RR] and halve the bound.
C3125 := HyperellipticCurve(x^5, R ! 1);  // 3125.a.3125.1: y^2 + y = x^5, CM by Q(zeta_5)
error if Sort(&cat RealRepresentationBound(C3125, 200)) ne ["CC", "CC"],
    "3125.a.3125.1: LMFDB factorsRR_geom is [CC, CC]";
