// Tests for the upperbounds Magma port.
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds.m Utils/Exit.m
// Assumes MAGMA_USER_SPEC points at CHIMP/CHIMP.spec, or attach manually.

AttachSpec("../endomorphisms/magma/spec");
AttachSpec("/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec");

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

// ----- EndomorphismAlgebraEtaBound -----
// Port of the eta/t narrowing step from Sage upper_bounds.py:55-78.
// Sage docstring example: F3, F7, F13 are Frobenius polys for a genus-2 AV
// that is geometrically isogenous to E^2 (so eta(A) = 4, t = 1).
// Coefficients are over Z (matches the L-polynomials emitted by LPolynomials).

ZZT<T> := PolynomialRing(Integers());
F3 := 1 - T^2 + 9*T^4;
F7 := 1 + 4*T^2 + 49*T^4;
F13 := 1 - 8*T^2 + 169*T^4;

ok, msg, eta_char0, t, eta_lower := EndomorphismAlgebraEtaBound([F3, F7, F13]);
assert ok;
assert eta_char0 eq 4;
assert t eq 1;
assert #eta_lower eq 3;
// Every entry of eta_lower should be a length-t endo factorization.
for endo in eta_lower do
    assert #endo eq t;
end for;

// User-provided eta_char0 hint: if the hint is correct (4), result should
// match. Internally eta := 2 * eta_char0 = 8.
ok2, _, eta_char0_2, t2, _ := EndomorphismAlgebraEtaBound([F3, F7, F13] : eta_char0 := 4);
assert ok2;
assert eta_char0_2 eq 4;
assert t2 eq 1;

// User-provided eta_char0 hint that is too small: no prime hits eta(A_p) = 2*eta(A),
// so we return success=false with eta_lower empty.
ok3, msg3, _, _, eta_lower3 := EndomorphismAlgebraEtaBound([F3, F7, F13] : eta_char0 := 1);
assert not ok3;
assert #eta_lower3 eq 0;

// ----- SubfieldsPolynomials -----
// Port of Sage subfields_polynomials. Returns polredabs'd defining polynomials
// of every subfield of NumberField(f), Q included. Already-canonical inputs
// are used so the test passes even without PARI/gp installed.

// Q(sqrt 5): subfields are Q and Q(sqrt 5).
sf := SubfieldsPolynomials(x^2 - 5);
assert #sf eq 2;
assert x in sf;
assert (x^2 - 5) in sf;

// Q(2^(1/4)): subfields are Q, Q(sqrt 2), Q(2^(1/4)).
sf := SubfieldsPolynomials(x^4 - 2);
assert #sf eq 3;
assert x in sf;
assert (x^2 - 2) in sf;
assert (x^4 - 2) in sf;

// ----- FieldIntersection -----
// Largest common subfield as a defining polynomial. Port of Sage
// field_intersection in endomorphisms/UpperBounds/utils.py.

// Isomorphic case: returns absolute polynomial of L (NOT polredabs'd; matches
// Sage). Q(sqrt 5) ~ Q(sqrt 5).
L := NumberField(x^2 - 5);
assert FieldIntersection(L, L) eq x^2 - 5;

// Linearly disjoint: Q(sqrt 2) and Q(sqrt 3) share only Q.
L := NumberField(x^2 - 2);
K := NumberField(x^2 - 3);
assert FieldIntersection(L, K) eq x;

// Nested: Q(sqrt 2) is a subfield of Q(2^(1/4)). Intersection is Q(sqrt 2).
L := NumberField(x^4 - 2);
K := NumberField(x^2 - 2);
assert FieldIntersection(L, K) eq x^2 - 2;

// ----- FieldIntersectionList -----
// Iterative intersection across a list of polynomials. Port of Sage
// field_intersection_list in endomorphisms/UpperBounds/utils.py.

// Singleton: intersection is the field itself.
assert FieldIntersectionList([x^2 - 5]) eq x^2 - 5;

// Linearly disjoint: only Q in common.
assert FieldIntersectionList([x^2 - 2, x^2 - 3]) eq x;

// Same field, two different defining polynomials.
assert FieldIntersectionList([x^2 - 2, x^2 - 8]) eq x^2 - 2;

// Any degree-1 polynomial collapses the intersection to Q.
assert FieldIntersectionList([x, x^2 - 3]) eq x;

// ----- FieldIntersectionMatrix -----
// Port of Sage field_intersection_matrix. For each column, compute the common
// subfields across rows (each row contributes the union of subfields of its
// polynomials). Output structure: <A_k, B_k> per column, where B_k is the
// sorted list of common-subfield defining polynomials and A_k is the maximal
// common subfield polynomial when uniquely determined (else the zero polynomial).

// Single-column shortcut: returns full subfield list of the intersection field.
M := [[x^2 - 5]];
result := FieldIntersectionMatrix(M);
assert #result eq 1;
assert result[1][1] eq x^2 - 5;
assert SequenceToSet(result[1][2]) eq {x, x^2 - 5};

// Single-column, linearly disjoint rows: intersection is Q.
M := [[x^2 - 2], [x^2 - 3]];
result := FieldIntersectionMatrix(M);
assert #result eq 1;
assert result[1][1] eq x;
assert result[1][2] eq [x];

// Multi-column, all rows agree.
M := [[x, x^2 - 2], [x, x^2 - 2]];
result := FieldIntersectionMatrix(M);
assert #result eq 2;
assert result[1][1] eq x;
assert result[1][2] eq [x];
assert result[2][1] eq x^2 - 2;
assert SequenceToSet(result[2][2]) eq {x, x^2 - 2};

// Multi-column, rows totally disagree: every column collapses to Q.
M := [[x^2 - 2, x^2 - 3], [x^2 - 5, x^2 - 7]];
result := FieldIntersectionMatrix(M);
assert #result eq 2;
assert result[1][1] eq x;
assert result[1][2] eq [x];
assert result[2][1] eq x;
assert result[2][2] eq [x];

// ----- EndomorphismAlgebraCenterBounds -----
// Port of Sage upper_bounds.py:82-125. Takes (eta, t, eta_lower) from EtaBound,
// checks multiset agreement across primes, and applies FieldIntersectionMatrix
// to bound the centers of each simple factor.

// Sage F3/F7/F13 example: drives the full output structure.
ZZT<T> := PolynomialRing(Integers());
F3 := 1 - T^2 + 9*T^4;
F7 := 1 + 4*T^2 + 49*T^4;
F13 := 1 - 8*T^2 + 169*T^4;
_, _, eta_c, t, eta_lower := EndomorphismAlgebraEtaBound([F3, F7, F13]);

ok, msg, output, total_dim := EndomorphismAlgebraCenterBounds(2 * eta_c, t, eta_lower);
assert ok;
assert total_dim eq 4;
assert #output eq 1;
ejnj, njdimAj, Lj, RRj := Explode(output[1]);
assert ejnj eq 2;
assert njdimAj eq 2;
assert RRj eq ["M_2(RR)"];
// The center is Q (intersection of three distinct imaginary quadratic fields).
assert Degree(Lj[2][#Lj[2]]) eq 1;

// Multiset disagreement case: build a synthetic eta_lower where the second
// prime has a different multiset shape. Should return success=false.
fake := [
    [<2, 4, T^2 + 1>],   // prime 1: one factor pair (2, 4)
    [<1, 2, T^2 + 1>, <1, 2, T^2 + 1>]  // prime 2: two factor pairs (1, 2)
];
ok2, msg2, _, _ := EndomorphismAlgebraCenterBounds(8, 1, fake);
assert not ok2;

print "Test-UpperBounds: all assertions passed.";
