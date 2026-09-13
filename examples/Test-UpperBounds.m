// Tests for the upperbounds Magma port.
// Run via: magma -b Utils/SetQuitOnErrortrue.m Test-UpperBounds.m Utils/Exit.m
// Assumes MAGMA_USER_SPEC points at CHIMP/CHIMP.spec, or attach manually.

AttachSpec("../endomorphisms/magma/spec");
AttachSpec("/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec");

R<x> := PolynomialRing(Rationals());

// ----- AlternatingSquare / SymmetricSquare CharacteristicPolynomial -----
// CMSV (arXiv:1705.09248) Section 7. Tensor(f,f) = alt^2 * f2, so alt has degree
// d*(d-1)/2. The literals below come from the elementary symmetric functions of
// f = (1 - a x)^2 (1 - b x)^2 with a + b = -1, a b = 2, not from the intrinsic.
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
    // Every kept prime must predict the enumerated point count.
    assert p + 1 + Coefficient(Lp, 1) eq #Points(ChangeRing(C, GF(p)));
end for;

// Cross-check: every (p, Lp) in the output must equal LPolynomial of the
// base-changed curve. Catches argument shuffling / convention errors. The
// algorithm is named so the check pins the pairing, not Magma's selection.
for pair in lpolys do
    p, Lp := Explode(pair);
    assert Lp eq LPolynomial(ChangeRing(C, GF(p)) : Al := "Naive");
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

// Magma's default algorithm returns a wrong (though formally valid) Weil
// polynomial at p = 5 here. That makes 5 the unique minimiser of
// dim End(Abar_p), so the bound names a quartic CM centre: the baseline read
// its real quadratic subfield and reported [RR, RR], and naming it now errors.
procedure test_lpolynomials_exact_at_small_primes()
    // LMFDB 8788.g.8788.1, geom_end_alg M_2(Q), factorsRR_geom [M_2(RR)].
    C8788 := HyperellipticCurve(R ! [-5, 10, -20, 19, -15, 6, -2], R ! [1, 1, 1]);

    // Point counts fix L_5 for genus 2, and are computed by enumeration, so
    // they are an oracle independent of any L-polynomial algorithm.
    C5 := ChangeRing(C8788, GF(5));
    assert #Points(C5) eq 6;
    assert #Points(BaseChange(C5, GF(25))) eq 28;
    lp5 := [pair[2] : pair in LPolynomials(C8788, 50) | pair[1] eq 5];
    assert #lp5 eq 1;
    error if lp5[1] ne 25*x^4 + x^2 + 1,
        Sprintf("8788.g.8788.1: #C(F_5) = 6 and #C(F_25) = 28 give L_5 = 25x^4 + x^2 + 1, got %o", lp5[1]);

    error if RealRepresentationBound(C8788, 50) ne [["M_2(RR)"]],
        Sprintf("8788.g.8788.1: LMFDB factorsRR_geom is [M_2(RR)], got %o",
                RealRepresentationBound(C8788, 50));
end procedure;

test_lpolynomials_exact_at_small_primes();

// ----- EndomorphismAlgebra over prime fields -----
// A Weil polynomial with constant term 1 has leading coefficient q^genus, so
// genus 1 needs the exact first root. Inputs are CLV (arXiv:1906.02803)
// Example 5.5: 11.a2 at p = 2, 3, where M(2) = Q(i) and M(3) = Q(sqrt(-11)).

ZZTe<Te> := PolynomialRing(Integers());
QTe<TQ> := PolynomialRing(Rationals());

// p = 2: supersingular, so End(Abar) is a quaternion algebra of dimension 4.
// M(2) = Q(i) pins Frobenius as -1 + i, whose fourth power is first rational,
// so endomorphisms appear over F_(2^4) and c = 4T + 1 with multiplicity 2.
procedure test_endomorphism_algebra_supersingular_prime_field()
    dim, fext, endo := EndomorphismAlgebra(1 + 2*Te + 2*Te^2);
    assert dim eq 4;
    assert fext eq 4;
    assert endo eq [<2, 2, 4*Te + 1>];
end procedure;

test_endomorphism_algebra_supersingular_prime_field();

// p = 3: ordinary, so End(Abar) is M(3) = Q(sqrt(-11)), dimension 2, already
// defined over F_3. c is compared as a field, since the paper fixes the field
// and not a defining polynomial.
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
// A_k bounds a centre, so every candidate must embed into it. Q(2^(1/4)) and
// Q(18^(1/4)) are non-isomorphic, incomparable, and share two distinct composita,
// so their common subfields have no greatest member. Needs gp for canonical names.
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

    // A third row equal to one quartic makes it the greatest candidate in every
    // row order. The pairwise fold keeps whichever it meets first, then collapses
    // to Q(sqrt 2), so it gets some of these four matrices wrong.
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

// No greatest candidate means no centre bound: CMSV Lemma 7.4.2 needs the true
// centre to embed in the reported field. Degree 8 is out of reach of a Frobenius
// factor below genus 4, so the entries are field-theoretic rather than Weil.
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
// Naming the CM field's real quadratic subfield instead reports [RR, RR] and
// halves the dimension. Expected values are LMFDB factorsRR_geom, independent
// of this implementation.
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

// ----- Zywina (arXiv:2009.07441) "A CM example": genus 4 -----
// y^2 = x^9 - 1. The endomorphisms of A_Qbar are defined over Q(zeta_9), so only
// p = 1 mod 9 sees them all; Shioda's A ~ B x E, with B simple of dimension 3
// carrying Z[zeta_9] and E elliptic, gives Q(zeta_9) x Q(zeta_3), dimension 8.
procedure test_upper_bound_zywina_cm_genus_4()
    C := HyperellipticCurve(x^9 - 1);
    eligible := [];
    for pair in LPolynomials(C, 100) do
        _, fieldext := EndomorphismAlgebra(pair[2]);
        if fieldext eq 1 then
            Append(~eligible, pair[1]);
        end if;
    end for;
    error if eligible ne [p : p in PrimesUpTo(100) | p mod 9 eq 1],
        Sprintf("x^9 - 1: only p = 1 mod 9 defines every endomorphism, got %o", eligible);

    ok, _, _, t, output, total_dim := EndomorphismAlgebraUpperBound(C, 100);
    assert ok;
    error if t ne 2, Sprintf("x^9 - 1: A ~ B x E has two factors, got %o", t);
    error if total_dim ne 8,
        Sprintf("x^9 - 1: Q(zeta_9) x Q(zeta_3) has dimension 8, got %o", total_dim);
    assert Sort([tup[2] : tup in output]) eq [1, 3];
    for tup in output do
        _, njdimAj, Lj, RRj := Explode(tup);
        n := njdimAj eq 3 select 9 else 3;
        error if not IsIsomorphic(NumberField(Lj[1]), CyclotomicField(n)),
            Sprintf("x^9 - 1: dimension %o factor has center Q(zeta_%o), got %o",
                    njdimAj, n, Lj[1]);
        error if RRj ne ["CC" : i in [1..njdimAj]],
            Sprintf("x^9 - 1: dimension %o factor is CM, got %o", njdimAj, RRj);
    end for;
end procedure;

test_upper_bound_zywina_cm_genus_4();

// ----- Zywina (arXiv:2009.07441) "Another example": genus 10 -----
// y^2 = x(x^20 + 7x^18 - 7x^2 - 1). End(A_Qbar) tensor Q is a definite quaternion
// algebra over Q: one factor, center Q, dimension 4, and tensor RR equal to HH.
// Zywina reads HH off explicit automorphisms of C, not off Frobenius data.
procedure test_upper_bound_zywina_quaternion_genus_10()
    C := HyperellipticCurve(x*(x^20 + 7*x^18 - 7*x^2 - 1));
    lpolys := LPolynomials(C, 50);
    // The set P of Section 1.5 meets [1, 50) in 17 and 41, and there
    // P_{A,p} = Q_p^2 with Q_p of degree 10 whose discriminant lies in -2 (Q^*)^2.
    for p in [17, 41] do
        fac := Factorization([pair[2] : pair in lpolys | pair[1] eq p][1]);
        error if [<Degree(f[1]), f[2]> : f in fac] ne [<10, 2>],
            Sprintf("p = %o: L_p should be an irreducible of degree 10, squared; got %o",
                    p, fac);
        error if not IsSquare(-Discriminant(fac[1][1]) / 2),
            Sprintf("p = %o: disc of the degree 10 factor should lie in -2 (Q^*)^2", p);
    end for;

    ok, _, _, t, output, total_dim := EndomorphismAlgebraUpperBound(
        [pair[2] : pair in lpolys]);
    assert ok;
    error if t ne 1, Sprintf("genus 10: A_Qbar is simple, got %o factors", t);
    error if total_dim ne 4,
        Sprintf("genus 10: the quaternion algebra has dimension 4, got %o", total_dim);
    assert #output eq 1;
    ejnj, njdimAj, Lj, RRj := Explode(output[1]);
    assert <ejnj, njdimAj> eq <2, 10>;
    error if Degree(Lj[1]) ne 1, Sprintf("genus 10: the center is Q, got %o", Lj[1]);
    // The truth is ["HH"]. Zywina states Frobenius polynomials cannot separate the
    // two, so the undecided pair is the sharpest sound answer this method can give.
    error if RRj ne ["M_2(RR) or HH"],
        Sprintf("genus 10: expected the type II/III ambiguity, got %o", RRj);
end procedure;

test_upper_bound_zywina_quaternion_genus_10();
