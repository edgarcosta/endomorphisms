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

print "Test-UpperBounds: all assertions passed.";
