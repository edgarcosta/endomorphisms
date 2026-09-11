// Tests for the upperbounds Magma port.
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds.m Utils/Exit.m
// Assumes MAGMA_USER_SPEC points at CHIMP/CHIMP.spec, or attach manually.

AttachSpec("../endomorphisms/magma/spec");
AttachSpec("/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec");

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

// ----- AlternatingSquare / SymmetricSquare CharacteristicPolynomial -----
// These are the exterior/symmetric square constructions of CMSV
// (arXiv:1705.09248) Section 7. The defining identity is
// TensorCharacteristicPolynomial(f, f) = alt^2 * f2 (f2 the char poly of the
// square), so alt has degree d*(d-1)/2 and the symmetric square Tensor/alt has
// degree d*(d+1)/2. A wrong degree here is exactly what the missing product in
// the degree assertion of AlternatingSquareCharacteristicPolynomial hid.

// f = (1 + x + 2x^2)^2 is the square of the Frobenius polynomial of an elliptic
// curve over F_2 with a_2 = -1. Writing 1 + x + 2x^2 = (1 - a x)(1 - b x), so
// a + b = -1 and a b = 2, the input has reciprocal roots a, a, b, b and its
// alternating square has reciprocal roots a^2, b^2 and a b with multiplicity 4.
// Since a^2 + b^2 = -3 and a^2 b^2 = 4, that product is
// (1 + 3x + 4x^2)(1 - 2x)^4, expanded below. The literal is therefore derived
// from the elementary symmetric functions of the input, not from the intrinsic.
procedure test_alternating_and_symmetric_square_degree_4()
    fsq := (1 + x + 2*x^2)^2;
    altsq := AlternatingSquareCharacteristicPolynomial(fsq);
    assert Degree(altsq) eq 6;
    assert altsq eq (1 + 3*x + 4*x^2) * (1 - 2*x)^4;
    assert altsq eq 64*x^6 - 80*x^5 + 16*x^4 + 8*x^3 + 4*x^2 - 5*x + 1;
    assert altsq^2 * PowerCharacteristicPolynomial(fsq, 2)
        eq TensorCharacteristicPolynomial(fsq, fsq);

    symsq := SymmetricSquareCharacteristicPolynomial(fsq);
    assert Degree(symsq) eq 10;
    assert altsq * symsq eq TensorCharacteristicPolynomial(fsq, fsq);
end procedure;

test_alternating_and_symmetric_square_degree_4();

// d = 2: the alternating square is the determinant, a linear factor.
procedure test_alternating_square_degree_2()
    assert Degree(AlternatingSquareCharacteristicPolynomial(1 + x + 2*x^2)) eq 1;
end procedure;

test_alternating_square_degree_2();

// d = 4: a genuine genus-2 Weil polynomial.
procedure test_alternating_square_genus_2_weil()
    assert Degree(AlternatingSquareCharacteristicPolynomial(1 - x^2 + 9*x^4)) eq 6;
end procedure;

test_alternating_square_genus_2_weil();

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

// ----- EndomorphismAlgebra over prime fields -----
// The leading coefficient of a Weil polynomial normalised to constant term 1
// is q^genus, so genus 1 over a prime field needs the exact 1st root; asking
// for a nontrivial perfect power rejects every genus-1 Weil polynomial over a
// prime field. The two genus-1 inputs are the Frobenius polynomials of CLV
// (arXiv:1906.02803) Example 5.5, elliptic curve 11.a2 at p = 2 and p = 3,
// where the paper records M(2) = Q(sqrt(-1)) and M(3) = Q(sqrt(-11)).

ZZTe<Te> := PolynomialRing(Integers());
QTe<TQ> := PolynomialRing(Rationals());

// p = 2: a_2 = -2 is even, so 11.a2 is supersingular at 2 and End(Abar) is a
// quaternion algebra, of dimension 4 over Q. The paper's M(2) = Q(sqrt(-1))
// pins Frobenius as alpha = -1 + i up to conjugacy; alpha^2 = -2i is irrational
// while alpha^4 = -4 is rational, so all endomorphisms first appear over
// F_(2^4) and det(1 - T Frob^4 | H^1) = (1 + 4T)^2. That is one factor
// <m, m * deg c, c> with m = 2 and c = 4T + 1, computed from the paper's field
// rather than read off the intrinsic.
procedure test_endomorphism_algebra_supersingular_prime_field()
    dim, fext, endo := EndomorphismAlgebra(1 + 2*Te + 2*Te^2);
    assert dim eq 4;
    assert fext eq 4;
    assert endo eq [<2, 2, 4*Te + 1>];
end procedure;

test_endomorphism_algebra_supersingular_prime_field();

// p = 3: a_3 = -1 is prime to 3, so 11.a2 is ordinary at 3, Abar stays simple,
// and End(Abar) is the imaginary quadratic field the paper records as
// M(3) = Q(sqrt(-11)), of dimension 2 over Q and already defined over F_3.
// The returned factor c is compared to Q(sqrt(-11)) as a field rather than to
// a particular defining polynomial: the paper fixes the field, not the model.
procedure test_endomorphism_algebra_ordinary_prime_field()
    dim, fext, endo := EndomorphismAlgebra(1 + Te + 3*Te^2);
    assert dim eq 2;
    assert fext eq 1;
    assert #endo eq 1;
    m, dimfactor, c := Explode(endo[1]);
    assert m eq 1;
    assert dimfactor eq 2;
    assert IsIsomorphic(NumberField(QTe!Reverse(c)), QuadraticField(-11));
end procedure;

test_endomorphism_algebra_ordinary_prime_field();

// Genus-2 control: the exact square root of 9 = 3^2 is still accepted.
procedure test_endomorphism_algebra_genus_2_square_leading_coefficient()
    dim, fext := EndomorphismAlgebra(1 - Te^2 + 9*Te^4);
    assert dim eq 8;
    assert fext eq 2;
end procedure;

test_endomorphism_algebra_genus_2_square_leading_coefficient();

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
// of every subfield of NumberField(f), Q included. The subfields themselves are
// the behavior being pinned; which polynomial names each one is decided by
// whichever polredabs is available (Polredabs returns a non-canonical fallback
// when PARI/gp is not on PATH), so subfields are compared up to isomorphism
// unless their polynomial is a polredabs fixed point.

// Q(sqrt 5): subfields are Q and Q(sqrt 5). Q(sqrt 5) is named by x^2 - x - 1
// with PARI/gp present and by x^2 - 5 without it.
sf := SubfieldsPolynomials(x^2 - 5);
assert #sf eq 2;
assert x in sf;
assert exists{p : p in sf | Degree(p) eq 2
                            and IsIsomorphic(NumberField(p), NumberField(x^2 - 5))};

// Q(2^(1/4)): subfields are Q, Q(sqrt 2), Q(2^(1/4)). All three polynomials are
// polredabs fixed points, so exact equality holds either way.
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

// Singleton: intersection is the field itself, Q(sqrt 5). The return value is
// polredabs'd, so pin the field rather than a defining polynomial.
inter5 := FieldIntersectionList([x^2 - 5]);
assert Degree(inter5) eq 2;
assert IsIsomorphic(NumberField(inter5), NumberField(x^2 - 5));

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
// A_1 and B_1 are polredabs'd, so Q(sqrt 5) is pinned up to isomorphism.
M := [[x^2 - 5]];
result := FieldIntersectionMatrix(M);
assert #result eq 1;
assert Degree(result[1][1]) eq 2;
assert IsIsomorphic(NumberField(result[1][1]), NumberField(x^2 - 5));
// B_1 is Q and Q(sqrt 5), sorted ascending by degree.
assert #result[1][2] eq 2;
assert result[1][2][1] eq x;
assert Degree(result[1][2][2]) eq 2;
assert IsIsomorphic(NumberField(result[1][2][2]), NumberField(x^2 - 5));

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

// ----- One entry per row: the whole candidate family, not a pairwise fold -----
// A_k is used as a center bound, so every member of B_k must embed into it, and
// B_k must therefore be the complete family of common subfields. A single-column
// matrix used to be answered by folding FieldIntersection pairwise, which keeps
// one common subfield of greatest degree per step and drops the alternatives, so
// it could name an A_k that other candidates do not embed into (issue
// endomorphisms-0t4).
//
// Witnesses: N = Q(2^(1/4), i, sqrt 3) is Galois over Q with group D_4 x C_2.
// Writing s for complex conjugation, r for the order-4 rotation in D_4 and c for
// the generator of C_2, the index-4 subgroups <s, c> and <s, r^2 c> are not
// conjugate, so the quartic fields Q(2^(1/4)) and Q(18^(1/4)), with
// 18^(1/4) = sqrt3 * 2^(1/4), are non-isomorphic; they are also incomparable,
// since Q(sqrt 2) is the only quadratic subfield of Q(2^(1/4)) and 18^(1/4)
// there would force sqrt 3 in as well. Their compositum is Q(2^(1/4), sqrt 3)
// for one embedding and Q(2^(1/4), sqrt -3) for the other, so those two octic
// fields share both quartics and nothing above them. That non-unique compositum
// is what makes a family of common subfields lack a greatest member, and degree 8
// is the smallest degree where it happens.
//
// Fields are compared through their Polredabs polynomials, so these scenarios
// need PARI/gp on PATH: with the non-canonical fallback no two rows agree on a
// name for the same field.
function d4c2_octic_pair()
    K := NumberField(x^4 - 2);
    return R ! DefiningPolynomial(AbsoluteField(ext<K | Polynomial([K | -3, 0, 1])>)),
           R ! DefiningPolynomial(AbsoluteField(ext<K | Polynomial([K | 3, 0, 1])>));
end function;

procedure test_field_intersection_matrix_incomparable_candidates()
    fX, fY := d4c2_octic_pair();
    Ak, Bk := Explode(FieldIntersectionMatrix([[fX], [fY]])[1]);
    // Q, Q(sqrt 2), Q(2^(1/4)) and Q(18^(1/4)) embed in both octic fields.
    error if [Degree(f) : f in Bk] ne [1, 2, 4, 4],
        Sprintf("expected candidates of degree [1, 2, 4, 4], got %o", Bk);
    error if not exists{f : f in Bk | IsIsomorphic(NumberField(f), NumberField(x^4 - 2))},
        Sprintf("Q(2^(1/4)) is missing from the candidates %o", Bk);
    error if not exists{f : f in Bk | IsIsomorphic(NumberField(f), NumberField(x^4 - 18))},
        Sprintf("Q(18^(1/4)) is missing from the candidates %o", Bk);
    // Neither quartic embeds in the other, so no candidate bounds the rest.
    error if not IsZero(Ak),
        Sprintf("expected the no-greatest-candidate marker 0, got %o (gp on PATH?)", Ak);

    // A third row equal to one of the two quartics resolves the family to that
    // quartic: it embeds in both octic fields and the other quartic does not
    // embed in it, so it is the greatest candidate whatever order the rows come
    // in. The pairwise fold keeps whichever quartic it meets first and then
    // collapses to Q(sqrt 2), which the center need not embed into, so it gets
    // some of these four matrices right and the others wrong.
    for J in [x^4 - 2, x^4 - 18] do
        for M in [[[fX], [fY], [J]], [[J], [fX], [fY]]] do
            Ak, Bk := Explode(FieldIntersectionMatrix(M)[1]);
            error if [Degree(f) : f in Bk] ne [1, 2, 4],
                Sprintf("expected candidates of degree [1, 2, 4] with third field %o, got %o",
                        J, Bk);
            error if Degree(Ak) ne 4 or not IsIsomorphic(NumberField(Ak), NumberField(J)),
                Sprintf("expected the greatest candidate %o, got %o", J, Ak);
        end for;
    end for;
end procedure;

test_field_intersection_matrix_incomparable_candidates();

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

// A factor whose candidates have no greatest member has no center bound at all:
// CMSV Lemma 7.4.2 requires the true center to embed in the reported field, and
// a candidate incomparable to it does not bound it. One factor per prime, of
// shape <m, m * deg(h)> = <1, 8>, is the one-entry-per-row case; the pairwise
// fold used to report whichever of the two quartics it met first as the center,
// with real representation [RR, RR, RR, RR] and dimension 4. The entries are the
// field-theoretic
// witnesses above rather than Weil polynomials, since degree 8 is out of reach
// of a Frobenius factor below genus 4.
procedure test_center_bounds_refuses_incomparable_candidates()
    fX, fY := d4c2_octic_pair();
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        2, 1, [[<1, 8, fX>], [<1, 8, fY>]]);
    error if ok,
        Sprintf("expected no center bound for incomparable candidates, got %o of dimension %o",
                output, total_dim);
end procedure;

test_center_bounds_refuses_incomparable_candidates();

// ----- EndomorphismAlgebraUpperBound + RealRepresentationBound (frob_list) -----
// Top-level wrappers. Sage docstring example must round-trip:
//   endomorphisms_upper_bound([[3,F3],[7,F7],[13,F13]]) ==
//     (True, 'We have...', 4, 1, [(2, 2, [T, [T]], ['M_2(RR)'])], 4)

ok, msg, eta_c2, t2, output2, total_dim2 := EndomorphismAlgebraUpperBound([F3, F7, F13]);
assert ok;
assert eta_c2 eq 4;
assert t2 eq 1;
assert total_dim2 eq 4;
assert #output2 eq 1;
ej, nj, Lj, RR := Explode(output2[1]);
assert ej eq 2;
assert nj eq 2;
assert RR eq ["M_2(RR)"];

// RealRepresentationBound returns just the RRj lists.
rr_only := RealRepresentationBound([F3, F7, F13]);
assert rr_only eq [["M_2(RR)"]];

// Failure passthrough: the eta_char0 = 1 path that EtaBound rejects.
ok_f, _, _, _, _, _ := EndomorphismAlgebraUpperBound([F3, F7, F13] : eta_char0 := 1);
assert not ok_f;

// ----- Curve-level overloads -----
// Compose LPolynomials with the frob_list form. Genus-2 example: a curve known
// (LMFDB 169.a.169.1) to be the geom-Jacobian of E^2 -> factorsRR_geom = [M_2(RR)].
// Curve: y^2 = x^5 + x^4 (this is the historical genus2_hyperelliptic[169] case).
C169 := HyperellipticCurve(x^5 + x^4, x^3 + x + 1);
ok_c, _, _, _, output_c, total_dim_c := EndomorphismAlgebraUpperBound(C169, 30);
assert ok_c;
assert total_dim_c eq 4;
assert #output_c eq 1;
assert output_c[1][4] eq ["M_2(RR)"];
assert RealRepresentationBound(C169, 30) eq [["M_2(RR)"]];

// No prime below B has good reduction: CLV (arXiv:1906.02803) Algorithm 5.1
// allows any B ge 1 and prescribes failure rather than an error, so both
// curve-level overloads must report failure instead of propagating the empty
// frob_list into the SeqEnum forms. B = 1 admits no primes at all.
procedure test_curve_level_overloads_without_good_reduction_prime()
    Cempty := HyperellipticCurve(x^3 + x^2, x^3 + 1);
    assert LPolynomials(Cempty, 1) eq [];
    ok, msg, eta, t, output, dim := EndomorphismAlgebraUpperBound(Cempty, 1);
    assert not ok;
    assert #msg gt 0;
    assert eta eq 0;
    assert t eq 0;
    assert output eq [];
    assert dim eq 0;
    assert RealRepresentationBound(Cempty, 1) eq [];
end procedure;

test_curve_level_overloads_without_good_reduction_prime();
// ----- Quartic CM stratum: the center bound is the CM field -----
// Pins which common subfield is selected as the center of a geometrically
// simple genus-2 Jacobian with quartic CM. The bound must name the quartic CM
// field, whose real representation is [CC, CC]; naming its real quadratic
// subfield instead reports [RR, RR] and halves the dimension bound. Expected
// values come from LMFDB factorsRR_geom for the labels below, an oracle
// independent of this implementation; the dimension 4 is the degree of the CM
// field, which for a simple CM abelian surface is the whole geometric
// endomorphism algebra. Background: CMSV arXiv:1705.09248 Section 7 and CLV
// arXiv:1906.02803.
procedure test_upper_bound_quartic_cm_center()
    // LMFDB 3125.a.3125.1: y^2 + y = x^5, with CM by Q(zeta_5).
    C := HyperellipticCurve(x^5, R ! 1);
    ok, _, _, t, output, total_dim := EndomorphismAlgebraUpperBound(C, 200);
    assert ok;
    assert t eq 1;
    rr := &cat [tup[4] : tup in output];
    error if rr ne ["CC", "CC"],
        Sprintf("3125.a.3125.1: LMFDB factorsRR_geom is [CC, CC], got %o", rr);
    error if total_dim ne 4,
        Sprintf("3125.a.3125.1: CM field has degree 4, got dimension %o", total_dim);

    // LMFDB 28561.a.371293.1: a second curve on the same stratum.
    C := HyperellipticCurve(R ! [-2, 3, 2, -2, -2], R ! [0, 0, 0, 1]);
    rr := RealRepresentationBound(C, 200);
    error if rr ne [["CC", "CC"]],
        Sprintf("28561.a.371293.1: LMFDB factorsRR_geom is [CC, CC], got %o", rr);
end procedure;

test_upper_bound_quartic_cm_center();

print "Test-UpperBounds: all assertions passed.";
