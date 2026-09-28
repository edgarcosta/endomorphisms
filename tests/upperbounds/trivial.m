AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

R<x> := PolynomialRing(Rationals());

procedure test_list_certificate_is_exact()
    P2 := x^4 + 2*x^3 + 3*x^2 + 4*x + 4;
    P7 := x^4 + x^3 - 2*x^2 + 7*x + 49;
    L7 := R ! Reverse(P7);

    ok, certificate := CanProveTrivialEndomorphismRing([P2, L7]);
    assert ok;
    assert certificate eq <false, 1, 2, P2, P7, P7>;
end procedure;

test_list_certificate_is_exact();

procedure test_insufficient_lists_return_false()
    P2 := x^4 + 2*x^3 + 3*x^2 + 4*x + 4;
    assert not CanProveTrivialEndomorphismRing([R |]);
    assert not CanProveTrivialEndomorphismRing([P2]);
    assert not CanProveTrivialEndomorphismRing([P2, P2]);
end procedure;

test_insufficient_lists_return_false();

procedure test_geometrically_trivial_curve_succeeds_in_both_modes()
    C := HyperellipticCurve(x^2 + x, x^3 + 1);
    assert CanProveTrivialEndomorphismRing(C, 8);
    assert CanProveTrivialEndomorphismRing(C, 8 : Geometric := true);
end procedure;

test_geometrically_trivial_curve_succeeds_in_both_modes();

procedure test_curve_skips_lower_genus_reductions()
    P2<X, Y, Z> := ProjectiveSpace(Rationals(), 2);
    C := Curve(P2, X^4 + 2*X^3*Y + 2*Y^4 + Y^3*Z +
        2*X^2*Z^2 + Y*Z^3 + 2*Z^4);
    assert CanProveTrivialEndomorphismRing(C, 8);
end procedure;

test_curve_skips_lower_genus_reductions();

procedure test_base_trivial_cyclotomic_curve_distinguishes_modes()
    C := HyperellipticCurve(x^5, R ! 1);
    assert CanProveTrivialEndomorphismRing(C, 12);
    assert not CanProveTrivialEndomorphismRing(C, 12 : Geometric := true);
end procedure;

test_base_trivial_cyclotomic_curve_distinguishes_modes();

procedure test_real_multiplication_curve_returns_false()
    C := HyperellipticCurve(-x^5, x^3 + x + 1);
    assert not CanProveTrivialEndomorphismRing(C, 32);
    assert not CanProveTrivialEndomorphismRing(C, 32 : Geometric := true);
end procedure;

test_real_multiplication_curve_returns_false();

procedure test_complex_multiplication_curve_returns_false()
    C := HyperellipticCurve(x^5 - x);
    assert not CanProveTrivialEndomorphismRing(C, 32);
    assert not CanProveTrivialEndomorphismRing(C, 32 : Geometric := true);
end procedure;

test_complex_multiplication_curve_returns_false();

procedure test_split_curve_returns_false()
    C := HyperellipticCurve(x^4 + x^2, x^3 + 1);
    assert not CanProveTrivialEndomorphismRing(C, 32);
    assert not CanProveTrivialEndomorphismRing(C, 32 : Geometric := true);
end procedure;

test_split_curve_returns_false();

procedure test_malformed_polynomial_raises_error()
    malformed := 1 + x + 4*x^4;
    caught := false;
    try
        ok := CanProveTrivialEndomorphismRing([malformed]);
    catch e
        caught := true;
    end try;
    assert caught;
end procedure;

test_malformed_polynomial_raises_error();

exit;
