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

intrinsic FieldIntersection(L::FldNum, K::FldNum) -> RngUPolElt
{Defining polynomial of one common subfield of L and K of greatest degree, or x
 for Q. Greatest degree is not greatest: with several incomparable maximal
 common subfields the one returned need not contain the others, so it bounds
 nothing. See also FieldIntersectionMatrix, which computes the whole family}
    if IsIsomorphic(L, K) then
        return DefiningPolynomial(AbsoluteField(L));
    end if;
    // Magma's Subfields(K) returns K-and-proper-subfields excluding Q.
    // Sort ascending by absolute degree so we can iterate large-to-small.
    Ksubs := [s[1] : s in Subfields(K)];
    Lsubs := [s[1] : s in Subfields(L)];
    Sort(~Ksubs, func<a, b | Degree(AbsoluteField(a)) - Degree(AbsoluteField(b))>);
    Sort(~Lsubs, func<a, b | Degree(AbsoluteField(a)) - Degree(AbsoluteField(b))>);
    Kupper := #Ksubs;
    for i := #Lsubs to 1 by -1 do
        sL := Lsubs[i];
        dL := Degree(AbsoluteField(sL));
        for j := Kupper to 1 by -1 do
            sK := Ksubs[j];
            dK := Degree(AbsoluteField(sK));
            if dK eq dL then
                if IsIsomorphic(sL, sK) then
                    return Polredabs(DefiningPolynomial(AbsoluteField(sL)));
                end if;
            elif dK lt dL then
                break;
            else
                Kupper -:= 1;
            end if;
        end for;
    end for;
    // No non-trivial common subfield: return the polynomial defining Q.
    return Parent(DefiningPolynomial(L)).1;
end intrinsic;

// Discriminant of the primitive integral model of f. For an irreducible integral
// polynomial, monic or not, that discriminant is index^2 times the discriminant
// of the number field it defines, so every prime ramified in the field divides it.
function integral_discriminant(f)
    d := LCM([Integers() | Denominator(Rationals() ! c) : c in Coefficients(f)]);
    return Discriminant(PrimitivePart(PolynomialRing(Integers()) ! (d * f)));
end function;

intrinsic FieldIntersectionMatrix(M::SeqEnum[SeqEnum[RngUPolElt]]) -> SeqEnum
{Per column k of M, the pair <A_k, B_k>: B_k lists every common subfield of that
 column across the rows, ascending by degree, and A_k defines the greatest member
 of B_k, or is zero when B_k has none. Unlike the Sage original this takes no
 single-column shortcut, which could drop maximal candidates from B_k}
    require #M gt 0: "matrix must have at least one row";
    require #M[1] gt 0: "matrix must have at least one column";
    R := Parent(M[1][1]);

    // One entry per row: a candidate above Q is ramified somewhere (Minkowski),
    // so it divides the gcd below. Trivial gcd settles the column with no
    // Subfields or Polredabs call. Several columns take unions, so it cannot apply.
    if forall{row : row in M | #row eq 1} then
        if Abs(Gcd([Integers() | integral_discriminant(row[1]) : row in M])) eq 1 then
            return [<R.1, [R.1]>];
        end if;
    end if;

    // Pick the row with the smallest discriminant-GCD as the working row.
    Di := 1;
    Dmin := 0;
    for i in [1..#M] do
        D := 0;
        for poly in M[i] do
            D := Gcd(D, Integers() ! Discriminant(poly));
        end for;
        if Dmin eq 0 or D lt Dmin then
            Dmin := D;
            Di := i;
        end if;
        if Dmin eq 1 then break; end if;
    end for;

    working_row := M[Di];
    subfields_union := [SequenceToSet([R | ]) : i in [1..#M]];
    subfields_cached := [false : i in [1..#M]];
    output := [];
    for k in [1..#working_row] do
        f := working_row[k];
        Lsub := SequenceToSet(SubfieldsPolynomials(f));
        for i in [1..#M] do
            if i ne Di then
                if not subfields_cached[i] then
                    U := SequenceToSet([R | ]);
                    for g in M[i] do
                        U := U join SequenceToSet(SubfieldsPolynomials(g));
                    end for;
                    subfields_union[i] := U;
                    subfields_cached[i] := true;
                end if;
                Lsub := Lsub meet subfields_union[i];
                if #Lsub eq 1 then break; end if;
            end if;
        end for;
        Lsub_list := Sort(SetToSequence(Lsub),
                          func<a, b | Degree(a) - Degree(b)>);
        Ak := R ! 0;
        if #Lsub_list eq 1 then
            Ak := R.1;
        elif Degree(Lsub_list[#Lsub_list]) ne Degree(Lsub_list[#Lsub_list - 1]) then
            // Unique largest subfield: verify it actually contains all the
            // others. Compared as sets, since candidates of equal degree are
            // ordered arbitrarily on both sides.
            if SequenceToSet(SubfieldsPolynomials(Lsub_list[#Lsub_list])) eq Lsub then
                Ak := Lsub_list[#Lsub_list];
            end if;
        end if;
        Append(~output, <Ak, Lsub_list>);
    end for;
    return output;
end intrinsic;

intrinsic FieldIntersectionList(polys::SeqEnum[RngUPolElt]) -> RngUPolElt
{Defining polynomial of a common subfield of NumberField(f) over all f in polys,
 or x for Q. The pairwise fold is greedy and order-dependent, so the result need
 not have greatest degree among the common subfields, nor contain them. See also
 FieldIntersectionMatrix, which computes the whole family}
    require #polys gt 0: "polys must not be empty";
    R := Parent(polys[1]);
    // Any degree-1 polynomial defines Q, so the intersection is forced to Q.
    if exists{f : f in polys | Degree(f) eq 1} then
        return R.1;
    end if;
    // Fast path: GCD of discriminants of the first ~1000 polynomials. If it is 1
    // the fields have no common subfield beyond Q, which is weaker than being
    // linearly disjoint.
    D := 0;
    for i in [1..Min(#polys, 1000)] do
        D := Gcd(D, Integers() ! Discriminant(polys[i]));
        if D eq 1 then
            return R.1;
        end if;
    end for;
    // Iterative intersection.
    L := NumberField(polys[1]);
    for f in polys do
        K := NumberField(f);
        inter := FieldIntersection(L, K);
        if Degree(inter) eq 1 then
            return R.1;
        end if;
        L := NumberField(inter);
    end for;
    return Polredabs(DefiningPolynomial(AbsoluteField(L)));
end intrinsic;

intrinsic SubfieldsPolynomials(f::RngUPolElt) -> SeqEnum[RngUPolElt]
{Return polredabs'd defining polynomials of every subfield of NumberField(f),
 including the rational subfield Q (defined by the polynomial x). Port of
 Sage subfields_polynomials in endomorphisms/UpperBounds/utils.py. The list is
 sorted ascending by degree, so its last entry defines NumberField(f) itself}
    R := Parent(f);
    if Degree(f) eq 1 then
        return [Polredabs(R.1)];
    end if;
    out := [Polredabs(R.1)];
    K := NumberField(f);
    for sub in Subfields(K) do
        Append(~out, Polredabs(DefiningPolynomial(AbsoluteField(sub[1]))));
    end for;
    // Subfields returns no particular order, and callers read the largest
    // subfield off the end of this list.
    return Sort(out, func<a, b | Degree(a) - Degree(b)>);
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

