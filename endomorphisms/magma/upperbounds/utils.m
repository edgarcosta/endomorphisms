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
    assert Degree(res) eq d(d-1) div 2;
    return res;
end intrinsic;

intrinsic SymmetricSquareCharacteristicPolynomial(f::RngUPolElt) -> RngUPolElt
    {return the characteristic polynomial of the induced linear transformation on the symmetric square}
    // we will do this by dividing by the factor corresponding to the alternating square
    d := Degree(f);
    fof := TensorCharacteristicPolynomial(f, f);
    f2 := PowerCharacteristicPolynomial(f, 2);
    g := fof/f2;
    assert IsOne(Denominator(g));
    bool, altsquaref := IsSquare(Parent(f)!g);
    assert bool;
    assert Degree(altsquaref) eq d*(d-1) div 2;
    res := Parent(f)! fof/ altsquaref;
    assert IsOne(Denominator(res));
    assert Degree(res) eq d*(d+1) div 2;
    return res;
end intrinsic;

intrinsic LPolynomials(C::Crv, B::RngIntElt) -> SeqEnum
{Return a sequence of <p, L_p> tuples for primes p < B where C has good reduction,
 with L_p the L-polynomial of the reduction at p (constant term 1, degree 2*Genus(C)).
 Port of Sage get_frob_list_HyperellipticCurve, curve-type-agnostic.}
    g := Genus(C);
    out := [];
    p := 2;
    while p lt B do
        try
            Cp := ChangeRing(C, GF(p));
            Lp := LPolynomial(Cp);
            if Degree(Lp) eq 2 * g then
                Append(~out, <p, Lp>);
            end if;
        catch e
            ;
        end try;
        p := NextPrime(p);
    end while;
    return out;
end intrinsic;

intrinsic RealRepresentationString(g::RngIntElt, K::FldNum, d::RngIntElt) -> SeqEnum[MonStgElt]
    {Encode End(A^n) tensor RR as a list of strings, one per simple component.
     g is the dimension of the simple factor, K its center, d the dimension of
     the endomorphism algebra over K. Port of the Sage RR_representation routine
     in endomorphisms/UpperBounds/utils.py. See Section 7 of the paper.}
    if HasComplexConjugate(K) and not IsTotallyReal(K) then
        n := Degree(K) div 2;
        KRR := "CC";
    else
        KRR := "RR";
        n := Degree(K);
        if d mod 2 eq 0 and g mod 2 eq 0 and g gt 3 then
            // Type II/III ambiguity: type III is excluded for g <= 3.
            KRR := "M_2(RR) or HH";
            d := d div 2;
        end if;
    end if;
    if d gt 1 then
        out := Sprintf("M_%o(%o)", d, KRR);
    else
        out := KRR;
    end if;
    return [out : i in [1..n]];
end intrinsic;

