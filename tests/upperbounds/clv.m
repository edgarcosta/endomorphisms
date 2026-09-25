// Regression tests for the CLV centre step and its disjunctive output token.
// Silent on success. Run from the repository root via: ./tests/run.sh clv

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

Qx<x> := PolynomialRing(Rationals());
possibility := "oneof{CC|RR+RR+RR+RR}";

function RadicalPolynomial(ds)
    K := NumberField(x^2 - ds[1]);
    for d in ds[2 .. #ds] do
        K := AbsoluteField(ext<K | Polynomial([K | -d, 0, 1])>);
    end for;
    return Qx ! DefiningPolynomial(K);
end function;

function d4c2_octic_pair()
    K := NumberField(x^4 - 2);
    return Qx ! DefiningPolynomial(
               AbsoluteField(ext<K | Polynomial([K | -3, 0, 1])>)),
           Qx ! DefiningPolynomial(
               AbsoluteField(ext<K | Polynomial([K | 3, 0, 1])>));
end function;

function ExtendField(f, d)
    K := NumberField(f);
    return Qx ! DefiningPolynomial(
        AbsoluteField(ext<K | Polynomial([K | -d, 0, 1])>));
end function;

function SameField(f, g)
    return Degree(f) eq Degree(g) and
           IsIsomorphic(NumberField(f), NumberField(g));
end function;

function CenterOutput(polynomial_rows)
    eta_lower := [
        [<1, 4, Reverse(f)> : f in row]
        : row in polynomial_rows
    ];
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        8, #polynomial_rows[1], eta_lower);
    assert ok;
    return output, total_dim;
end function;

// Dropping an inclusion-maximal survivor because another survivor has larger
// degree loses the true imaginary-quadratic possibility in this fixture.
procedure test_astra_counterexample_keeps_all_literal_maxima()
    mq := RadicalPolynomial([2, 3, -1]);
    nq := RadicalPolynomial([2, 3, -13]);
    ap := RadicalPolynomial([5, 7, -1]);
    bp := RadicalPolynomial([2, 3, -11]);
    eta_lower := [
        [<1, 8, Reverse(mq)>, <1, 8, Reverse(nq)>],
        [<1, 8, Reverse(ap)>, <1, 8, Reverse(bp)>]
    ];

    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        16, 2, eta_lower);
    assert ok;
    assert total_dim eq 8;
    assert #output eq 2;
    assert output[1][3][1] eq 4;
    assert Sort([Degree(f) : f in output[1][3][2]]) eq [2, 4];
    assert output[1][4] eq [possibility];
    assert RealRepresentationEmbeds(["CC"], output[1][4]);
    assert RealRepresentationEmbeds(
        ["RR", "RR", "RR", "RR"], output[1][4]);
end procedure;

// Enumerating only irreducible factors misses the normic divisor x^4 + 2:
// over Q(2^(1/4)) it is a product of two non-normic quadratics.
procedure test_reducible_normic_divisor_in_non_disjoint_pair()
    eta_lower := [
        [<1, 4, Reverse(x^4 - 2)>],
        [<1, 4, Reverse(x^4 + 2)>]
    ];
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        8, 1, eta_lower);
    assert ok;
    assert total_dim eq 1;
    assert #output eq 1;
    assert output[1][3] eq <1, [x]>;
    assert output[1][4] eq ["RR"];
end procedure;

// Isomorphism reduction or a unique-greatest-field requirement discards one
// of these two distinct literal quartic survivors in the octic anchor.
procedure test_octic_pair_preserves_literal_maxima()
    fX, fY := d4c2_octic_pair();
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        2, 1, [[<1, 8, fX>], [<1, 8, fY>]]);
    assert ok;
    assert total_dim eq 4;
    assert #output eq 1;
    assert output[1][3][1] eq 4;
    assert #output[1][3][2] eq 2;
    polynomials := output[1][3][2];
    assert exists{f : f in polynomials | SameField(f, x^4 - 2)};
    assert exists{f : f in polynomials | SameField(f, x^4 - 18)};
    assert output[1][4] eq ["RR", "RR", "RR", "RR"];
end procedure;

// A third degree-eight field selects either quartic candidate in either
// anchor order. Pairwise intersection cannot preserve all four cases.
procedure test_octic_third_prime_selects_each_literal_quartic()
    fX, fY := d4c2_octic_pair();
    for J in [x^4 - 2, x^4 - 18] do
        fZ := ExtendField(J, 5);
        for polynomials in [[fX, fY, fZ], [fZ, fX, fY]] do
            ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
                2, 1, [[<1, 8, Reverse(f)>] : f in polynomials]);
            assert ok and total_dim eq 4 and #output eq 1 and
                   output[1][3][1] eq 4 and #output[1][3][2] eq 1 and
                   SameField(output[1][3][2][1], J);
        end for;
    end for;
end procedure;

// The S3 splitting field contains three conjugate cubic subfields. They are
// isomorphic but distinct literal subspaces and must survive every prime.
procedure test_conjugate_literal_subfields_survive_another_prime()
    f := x^3 - 4*x + 1;
    K := NumberField(f);
    m := Qx ! DefiningPolynomial(SplittingField(f));
    n2 := Qx ! DefiningPolynomial(AbsoluteField(
        ext<K | Polynomial([K | -2, 0, 1])>));
    n3 := Qx ! DefiningPolynomial(AbsoluteField(
        ext<K | Polynomial([K | -3, 0, 1])>));

    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        6, 1, [[<1, 6, Reverse(h)>] : h in [m, n2]]);
    assert ok and total_dim eq 3 and output[1][3][1] eq 3 and
           #output[1][3][2] eq 3;

    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        6, 1, [[<1, 6, Reverse(h)>] : h in [m, n2, n3]]);
    assert ok and total_dim eq 3 and output[1][3][1] eq 3 and
           #output[1][3][2] eq 3;
end procedure;

// Public center-bound cases replacing the relevant deleted helper contracts.
procedure test_center_selection_behaviours()
    output, _ := CenterOutput([[x^2 - 2], [x^2 - 2], [x^2 - 2]]);
    assert #output eq 1 and output[1][3][1] eq 2 and
           #output[1][3][2] eq 1 and
           SameField(output[1][3][2][1], x^2 - 2);

    output, _ := CenterOutput([[x^4 - 2], [x^2 - 2]]);
    assert #output eq 1 and output[1][3][1] eq 2 and
           #output[1][3][2] eq 1 and
           SameField(output[1][3][2][1], x^2 - 2);

    output, _ := CenterOutput([[x^2 - 5]]);
    assert #output eq 1 and output[1][3][1] eq 2 and
           #output[1][3][2] eq 1 and
           SameField(output[1][3][2][1], x^2 - 5);

    output, _ := CenterOutput([
        [x^2 - 2, x^2 - 3],
        [x^2 - 2, x^2 - 3]
    ]);
    assert #output eq 2 and
           SameField(output[1][3][2][1], x^2 - 2) and
           SameField(output[2][3][2][1], x^2 - 3);

    output, _ := CenterOutput([
        [x^2 - 2, x^2 - 3],
        [x^2 - 5, x^2 - 7]
    ]);
    assert #output eq 2 and
           forall{entry : entry in output | entry[3] eq <1, [x]>};
end procedure;

// CLV Example 5.6: retaining only the anchor leaves degree 2, while adding
// p = 37 cuts the centre degree to 1 and exposes the genus-4 ambiguity.
procedure test_clv_genus_four_example_degree_and_dimension()
    gs := [1 - 2*x + 19*x^2, 1 + 7*x + 37*x^2];

    ok, _, eta, t, eta_lower := EndomorphismAlgebraEtaBound([gs[1]^4]);
    assert ok and eta eq 16 and t eq 1;
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        2 * eta, t, eta_lower);
    assert ok;
    assert output[1][3][1] eq 2;
    assert total_dim eq 32;

    ok, _, eta, t, eta_lower := EndomorphismAlgebraEtaBound(
        [g^4 : g in gs]);
    assert ok and eta eq 16 and t eq 1;
    ok, _, output, total_dim := EndomorphismAlgebraCenterBounds(
        2 * eta, t, eta_lower);
    assert ok;
    assert output[1][3][1] eq 1;
    assert total_dim eq 16;
    assert output[1][4] eq ["oneof{M_4(RR)|M_2(HH)}"];
    assert RealRepresentationEmbeds(["M_2(HH)"], output[1][4]);
end procedure;

failures := [];
try
    test_astra_counterexample_keeps_all_literal_maxima();
catch e
    Append(~failures, <"astra", Sprint(e`Object)>);
end try;
try
    test_reducible_normic_divisor_in_non_disjoint_pair();
catch e
    Append(~failures, <"reducible divisor", Sprint(e`Object)>);
end try;
try
    test_octic_pair_preserves_literal_maxima();
catch e
    Append(~failures, <"octic pair", Sprint(e`Object)>);
end try;
try
    test_octic_third_prime_selects_each_literal_quartic();
catch e
    Append(~failures, <"octic third prime", Sprint(e`Object)>);
end try;
try
    test_conjugate_literal_subfields_survive_another_prime();
catch e
    Append(~failures, <"conjugate literal subfields", Sprint(e`Object)>);
end try;
try
    test_center_selection_behaviours();
catch e
    Append(~failures, <"center selection behaviours", Sprint(e`Object)>);
end try;
try
    test_clv_genus_four_example_degree_and_dimension();
catch e
    Append(~failures, <"CLV genus 4", Sprint(e`Object)>);
end try;
error if #failures ne 0, Sprintf("CLV regressions failed: %o", failures);
