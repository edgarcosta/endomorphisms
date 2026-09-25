declare verbose UpperBounds, 1;

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
    assert Degree(res) eq d*(d-1) div 2;
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
 L_p comes from exact point counting where affordable, and is kept only when it
 predicts #C(F_p). Port of Sage get_frob_list_HyperellipticCurve}
    return LPolynomials(C, 2, B);
end intrinsic;

intrinsic LPolynomials(C::Crv, Blow::RngIntElt, Bhigh::RngIntElt) -> SeqEnum
{Range form of LPolynomials: the <p, L_p> tuples for the good primes p with
 Blow <= p < Bhigh. A caller climbing a ladder of bounds can accumulate the
 tuples rung by rung and count points at each prime exactly once}
    g := Genus(C);
    out := [];
    p := Max(Blow, 2);
    if not IsPrime(p) then
        p := NextPrime(p);
    end if;
    while p lt Bhigh do
        try
            Cp := ChangeRing(C, GF(p));
            // Over small prime fields Magma's default returns a wrong but
            // formally valid Weil polynomial for some genus-2 curves. Count
            // points exactly while affordable (cutoff from UseZetaMethod);
            // LPolynomial caches on Cp and ignores Al after the first call.
            if p^g le 10^6 then
                Lp := LPolynomial(Cp : Al := "Naive");
            else
                Lp := LPolynomial(Cp);
            end if;
            if Degree(Lp) eq 2 * g then
                // Enumerating on an untouched reduction gives a trace oracle
                // owing nothing to either algorithm or to Lp's cache. Dropping
                // a prime that fails it only weakens the bound, where keeping a
                // wrong L_p would make it unsound. It checks the trace alone.
                n1 := #Points(ChangeRing(C, GF(p)));
                if p + 1 + Coefficient(Lp, 1) eq n1 then
                    Append(~out, <p, Lp>);
                else
                    vprintf UpperBounds: "p = %o dropped: L_p gives %o points, %o counted\n",
                        p, p + 1 + Coefficient(Lp, 1), n1;
                end if;
            end if;
        catch e
            ;
        end try;
        p := NextPrime(p);
    end while;
    return out;
end intrinsic;

intrinsic RealRepresentationString(g::RngIntElt, K::FldNum, d::RngIntElt) -> SeqEnum[MonStgElt]
    {The strings encoding End(A^n) tensor RR, one per simple component. Here g is
     the dimension n_j * dim(A_j) of the POWER, not of the simple factor, K is
     the center, and d is the DEGREE e_j * n_j of the central simple algebra over
     K, whose dimension over K is d^2. See arXiv:1705.09248, Section 7}
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

