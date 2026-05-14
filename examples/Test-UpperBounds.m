// Tests for the upperbounds Magma port.
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds.m Utils/Exit.m
// Assumes MAGMA_USER_SPEC points at CHIMP/CHIMP.spec, or attach manually.

AttachSpec("../endomorphisms/magma/spec");

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

// ----- RealRepresentationString -----
// Port of Sage RR_representation in endomorphisms/UpperBounds/utils.py:159-181.
// Encodes End(A^n) tensor RR as a list of strings, one per simple component.

// CM quadratic: K = Q(i), CM => "CC", n = deg(K)/2 = 1, d = 1.
K := NumberField(x^2 + 1);
assert RealRepresentationString(1, K, 1) eq ["CC"];

// Real quadratic, g = 2, d = 1: not CM, n = deg(K) = 2, KRR = "RR", d = 1 => "RR".
K := NumberField(x^2 - 2);
assert RealRepresentationString(2, K, 1) eq ["RR", "RR"];

// Real quadratic, g = 2, d = 2: not CM. d%2 = 0, g%2 = 0, but g <= 3 so
// the type II/III ambiguity rule does not fire. d > 1 => "M_2(RR)".
K := NumberField(x^2 - 2);
assert RealRepresentationString(2, K, 2) eq ["M_2(RR)", "M_2(RR)"];

// Real quadratic, g = 4, d = 2: not CM, d%2 = 0, g%2 = 0, g > 3 =>
// type II/III ambiguity: "M_2(RR) or HH", d collapses to 1.
K := NumberField(x^2 - 2);
assert RealRepresentationString(4, K, 2) eq ["M_2(RR) or HH", "M_2(RR) or HH"];

// ----- LPolynomials -----
// Port of Sage get_frob_list_HyperellipticCurve, curve-type-agnostic.
// Returns [<p, L_p> : p < B with good reduction].

// Hyperelliptic genus 2: y^2 = x^5 + x + 1
C := HyperellipticCurve(x^5 + x + 1);
g := Genus(C);
lpolys := LPolynomials(C, 20);
assert #lpolys gt 0;
for pair in lpolys do
    p, Lp := Explode(pair);
    assert IsPrime(p);
    assert p lt 20;
    assert Degree(Lp) eq 2 * g;
    assert Coefficient(Lp, 0) eq 1;
    assert Coefficient(Lp, 2 * g) eq p^g;
end for;

// Cross-check: every (p, Lp) in the output must equal LPolynomial of the
// base-changed curve. Catches argument shuffling / convention errors.
for pair in lpolys do
    p, Lp := Explode(pair);
    assert Lp eq LPolynomial(ChangeRing(C, GF(p)));
end for;

// Plane quartic genus 3: a smooth plane quartic.
P2<X, Y, Z> := ProjectiveSpace(Rationals(), 2);
Cpl := Curve(P2, X^4 + Y^4 + Z^4 + X*Y*Z*(X+Y+Z));
gpl := Genus(Cpl);
assert gpl eq 3;
lpolys_pl := LPolynomials(Cpl, 12);
assert #lpolys_pl gt 0;
for pair in lpolys_pl do
    p, Lp := Explode(pair);
    assert Degree(Lp) eq 2 * gpl;
    assert Coefficient(Lp, 0) eq 1;
end for;

print "Test-UpperBounds: all assertions passed.";
