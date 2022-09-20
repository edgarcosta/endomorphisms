intrinsic TensorCharacteristicPolynomial(f::RngUPolElt, g::RngUPolElt) -> RngUPolElt
    {given the characteristic polynomials of two linear transformations,
    return the characteristic polynomial of the induced linear transformation on the tensor product}
    require Parent(g) eq Parent(f) : "both arguments should have the same parent";
    R := Parent(g);
    k := BaseRing(R);
    _<x, y> := PolynomialRing(k, 2);
    A := Evaluate(f, y);
    B := Homogenization( Evaluate(g, x), y);
    return R!Coefficients(Resultant(A, B, y), x);
end intrinsic;

// Reverse(WeilPolynomialOverFieldExtension(Reverse(f), 2))
// if f is monic, then this matches WeilPolynomialOverFieldExtension(f, k)
// otherwise Reverse(WeilPolynomialOverFieldExtension(Reverse(f), k))
intrinsic PowerCharacteristicPolynomial(f::RngUPolElt, k::RngIntElt) -> RngUPolElt
    {return the characteristic polynomial of the kth power of the linear transformation}
    R := Parent(f);
    S<t, u> := PolynomialRing(BaseRing(R), 2);
    return R!Coefficients(Resultant(Evaluate(f, u), u^k - t, u), t);
end intrinsic;

intrinsic AlternatingSquareCharacteristicPolynomial(f::RngUPolElt) -> RngUPolElt
    {return the characteristic polynomial of the induced linear transformation on the alternating square}
    d := Degree(f);
    g := TensorCharacteristicPolynomial(f, f)/PowerCharacteristicPolynomial(f, 2);
    assert IsOne(Denominator(g));
    bool, res := IsSquare(Parent(f)!g);
    assert bool;
    assert Degree(res) eq d(d-1)/2;
    return res;
end intrinsic;

intrinsic SymmetricSquareCharacteristicPolynomial(f::RngUPolElt) -> RngUPolElt
    {return the characteristic polynomial of the induced linear transformation on the symmetric square}
    // we will do this by dividing by the factor corresponding to the alternating square
    fof := TensorCharacteristicPolynomial(f, f);
    f2 := PowerCharacteristicPolynomial(f, 2);
    g := fof/g;
    assert IsOne(Denominator(g));
    bool, altsquaref := IsSquare(Parent(f)!g);
    assert bool;
    assert Degree(altsquaref) eq d(d-1)/2;
    res := Parent(f)! fof/ altsquaref;
    assert IsOne(Denominator(res));
    assert Degree(res) eq d(d+1)/2;
    return res;
end intrinsic;

