/***
 *  Relative number field functionality and attributes
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


declare attributes FldNum : iota;
declare attributes FldRat : iota;


intrinsic IsQQ(F::Fld) -> BoolElt
{Returns whether or not the field F equals QQ.}

return Type(F) eq FldRat;

end intrinsic;


intrinsic HasBaseQQ(F::Fld) -> BoolElt
{Returns whether or not the field F has QQ as a base ring.}

return Type(BaseRing(F)) eq FldRat;

end intrinsic;


//intrinsic GiveName(~K::Fld, F::Fld, str::MonStgElt)
//{Give the name str to a generator of K over F if such a generator is indeed
//present.}
//
//if IsRelativeExtension(K, F) then
//    AssignNames(~K, [ str ]);
//end if;
//
//end intrinsic;


//intrinsic IsRelativeExtension(K::Fld, F::Fld) -> Fld
//{Returns whether K is an extension of F. Equalities are not considered as
//extensions.}
//
///* Over QQ, the field QQ is not an extension */
//testrat := IsQQ(F) and IsQQ(K);
///* Over general fields, F is not an extension, and neither is a field with a different base */
//testnonrat := (not IsQQ(F)) and (BaseRing(K) ne F);
///* Everything else is */
//return not (testrat or testnonrat);
//
//end intrinsic;


//intrinsic MakeRelative(K::Fld, F::Fld) -> Fld
//{Returns K as a relative field over F, taking a linear extension if K equals F.}
//
//if not IsRelativeExtension(K, F) then
//    R<x> := PolynomialRing(F);
//    K := NumberField(x - 1: DoLinearExtension := true);
//end if;
//return K;
//
//end intrinsic;


intrinsic FieldDescription(K::Fld) -> SeqEnum
{Returns a list describing the field K.}

/* Make relative to deal with the case of a linear extension */
F := BaseRing(K);
/* Use fewer brackets when the base field is QQ */
if IsQQ(F) then
    K_seq := [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
    K_seq := [ [ Rationals() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1)) ];
end if;
return K_seq;

end intrinsic;


intrinsic ElementDescription(r::.) -> .
{Returns a list describing the field element r.}

K := Parent(r); F := BaseRing(K);
/* Use fewer brackets when possible in the next two cases */
if IsQQ(K) then
    return r;
elif (IsQQ(F) and not IsQQ(K)) then
    return [ Rationals() ! c : c in Eltseq(r) ];
else
    return [ [ Rationals() ! c : c in Eltseq(d) ] : d in Eltseq(r) ];
end if;

end intrinsic;


intrinsic SetInfinitePlaceDownwardsRecursively(K::Fld, iota::.)
{Sets the infinite place of K to equal iota. This change is passed on
compatibly to the tower of fields to which K belongs.}

K`iota := iota;
if IsQQ(K) then
    return;
end if;
F := BaseRing(K);
/* The next clause is clumsy, but needed because of problems when encountering
 * the rational field when it is left out */
if IsQQ(F) then
    F`iota := InfinitePlaces(F)[1];
    return;
end if;
for iotaF in InfinitePlaces(F) do
    if Extends(iota, iotaF) then
        SetInfinitePlaceDownwardsRecursively(F, iotaF);
    end if;
end for;

end intrinsic;


intrinsic SetInfinitePlaceUpwards(K::Fld)
{Assigns an infinite place to K that is compatible with the infinite place of
the base field of K.}

F := BaseRing(K);
if not assigned F`iota or IsQQ(F) then
    SetInfinitePlaceDownwardsRecursively(K, InfinitePlaces(K)[1]);
else
    for iotaK in InfinitePlaces(K) do
        if Extends(iotaK, F`iota) then
            K`iota := iotaK;
        end if;
    end for;
end if;

end intrinsic;


intrinsic SetInfinitePlaceDownwards(L::Fld, K::Fld)
{Assigns an infinite place to K that is compatible with the infinite place of
the extension L of K. Both L and K should have the same base ring.}

if not assigned L`iota then
    return;
end if;
if IsQQ(K) then
    SetInfinitePlaceDownwardsRecursively(K, InfinitePlaces(K)[1]); return;
else
    /* TODO: This invocation is a bit strange, see if it could give rise to
     * bugs */
    test := IsSubfield(K, L);
    for iotaK in InfinitePlaces(K) do
        if Extends(L`iota, iotaK) then
            K`iota := iotaK; return;
        end if;
    end for;
end if;

end intrinsic;


intrinsic RationalsExtra() -> FldNum
{Returns the number field defined by f along with an infinite place.}

K := Rationals();
SetInfinitePlaceUpwards(K);
return K;

end intrinsic;


intrinsic NumberFieldExtra(f::RngUPolElt) -> FldNum
{Returns the number field defined by f along with an infinite place that
extends the place of the base field, if it has any.}

K := NumberField(f : DoLinearExtension := true);
SetInfinitePlaceUpwards(K);
return K;

end intrinsic;


intrinsic EmbedAtInfinitePlace(f::RngUPolElt, RCC::RngUPol) -> RngUPolElt
{Returns the polynomial f considered as a complex polynomial to precision
prec.}

R := Parent(f); d := #GeneratorsSequence(R);
F := BaseRing(R);
/* TODO: The field gets an infinite place if it did not have any, a slight
 * containment breach where side effects are concerned. One can be picky about
 * this in the future. */
if not assigned F`iota then
    SetInfinitePlaceDownwardsRecursively(F, InfinitePlaces(F)[1]);
end if;

if IsZero(f) then
    return RCC ! 0;
else
    F := BaseRing(f);
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ Evaluate(MonomialCoefficient(f, mon), F`iota : Precision := prec) * RCC.1^Degree(mon) : mon in mons ];
end if;

end intrinsic;


intrinsic EmbedAtInfinitePlace(f::RngMPolElt, RCC::RngMPol) -> RngMPolElt
{Returns the polynomial f considered as a complex polynomial to precision
prec.}

R := Parent(f); d := #GeneratorsSequence(R);
F := BaseRing(R);
/* TODO: The field gets an infinite place if it did not have any, a slight
 * containment breach where side effects are concerned. One can be picky about
 * this in the future. */
if not assigned F`iota then
    SetInfinitePlaceDownwardsRecursively(F, InfinitePlaces(F)[1]);
end if;

if IsZero(f) then
    return RCC ! 0;
else
    F := BaseRing(Parent(f));
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ Evaluate(MonomialCoefficient(f, mon), F`iota : Precision := prec) * Monomial(RCC, Exponents(mon)) : mon in mons ];
end if;

end intrinsic;


intrinsic EmbedAtInfinitePlace(fs::SeqEnum, prec::RngIntElt) -> SeqEnum
{Returns the list of polynomials fs considered as complex polynomials to
precision prec.}

R := Parent(fs[1]); d := #GeneratorsSequence(R);
F := BaseRing(R);
/* TODO: The field gets an infinite place if it did not have any, a slight
 * containment breach where side effects are concerned. One can be picky about
 * this in the future. */
if not assigned F`iota then
    SetInfinitePlaceDownwardsRecursively(F, InfinitePlaces(F)[1]);
end if;

CC<I> := ComplexFieldExtra(prec);
if d eq 1 then
    RCC := PolynomialRing(CC);
else
    RCC := PolynomialRing(CC, d);
end if;
fsCC := [ EmbedAtInfinitePlace(f, RCC) : f in fs ];
return fsCC;

end intrinsic;


intrinsic ClearDenominator(K::Fld) -> Fld
{Returns a field isomorphic to K over its base ring such that the defining
polynomial of the new field is integral. (Works over a larger class of base
fields than Polredabs.)}

F := BaseRing(K); d := Degree(K);
if d eq 1 then
    return K;
end if;
r := K.1; coeffs := Coefficients(MinimalPolynomial(r));
dens := Reverse([ Denominator(coeff) : coeff in coeffs ]); dens := dens[2..#dens];
primes := &join[ Set([ tup[1] : tup in Factorization(den) ]) : den in dens | den ne 0 ];
if #primes eq 0 then
    common_den := 1;
else
    common_den := &*[ p^Maximum([ Ceiling(Valuation(dens[k], p)/k) : k in [1..#dens] | dens[k] ne 0 ]) : p in primes ];
end if;

R := PolynomialRing(F);
f := MinimalPolynomial(common_den*r, F);
K := NumberField(R ! f);
return K;

end intrinsic;


intrinsic ImproveField(K::Fld) -> Fld
{Tries to realize to field K in as simple a way as possible.}

/*
K := ClearDenominator(K);
if HasBaseQQ(K) and not IsQQ(K) then
    return Polredabs(K);
end if;
*/

if IsQQ(K) then
    return K;
end if;

K := ClearDenominator(K);
if HasBaseQQ(K) then
    return Polredabs(K);
end if;
L := AbsoluteField(K); L := Polredabs(L);
test, h := IsSubfield(BaseRing(K), L);
return RelativeField(BaseRing(K), L);

end intrinsic;


intrinsic ExtendRelativeNumberField(f::RngUPolElt) -> Fld
{Given a polynomial f over a relative field K | F, returns an extension of F that contains both K and a root of f.}

K := BaseRing(f); F := BaseRing(K);
if Degree(f) eq 1 then
    return K;
end if;

vprintf EndoFind : "Extending %o by root of %o over %o...\n", K, f, F;
K := NumberField(f);
test, h := IsSubfield(F, K);
K := RelativeField(F, K);
/* TODO: The improvement costs time in general */
K := ImproveField(K);
vprintf EndoFind : "done\n";
return K;

end intrinsic;


intrinsic ExtendRelativeNumberFieldExtra(f::RngUPolElt) -> Fld
{Given a polynomial f over a relative field K | F, returns an extension of F
that contains both K and a root of f.}

K := ExtendRelativeNumberField(f);
SetInfinitePlaceUpwards(K);
return K;

end intrinsic;


intrinsic ExtendRelativeSplittingField(f::RngUPolElt) -> Fld
{Given a polynomial f over a relative field K | F, returns an extension of F
that contains both K and the splitting field of f.}

K := BaseRing(f); F := BaseRing(K);
vprintf EndoFind : "Extending %o by splitting %o over %o\n...", K, f, F;
if IsQQ(F) then
    if IsQQ(K) then
        K := ImproveField(SplittingField(Polredabs(f*LCM([ Integers() ! c : c in Coefficients(f) ]))));
    else
        K := ImproveField(SplittingField(f));
    end if;
    vprintf EndoFind : "done\n";
    return K;
end if;

while true do
    factors := [ tup[1] : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
    if #factors eq 0 then
        vprintf EndoFind : "done\n";
        return K;
    end if;
    K := ExtendRelativeNumberField(factors[1]);
end while;

end intrinsic;


function ComparePolynomials(f1, f2);
// Input:   Two subfields, fields, or polynomials.
// Output:  A comparison function: field with smaller degrees are smaller.

if Degree(f1) lt Degree(f2) then
    return -1;
elif Degree(f1) eq Degree(f2) then
    return 0;
else
    return 1;
end if;

end function;


intrinsic RelativeSplittingField(fs::SeqEnum : AssumeIrr := false) -> FldNum
{Returns a splitting field of the polynomials in fs over their common base
ring.}

F := BaseRing(fs[1]); R<x> := PolynomialRing(F);
K := NumberField(x - 1: DoLinearExtension := true);
//F := BaseRing(fs[1]); K := F;
fs := Reverse(Sort(fs, ComparePolynomials));
if not AssumeIrr then
    for f in fs do
        Fac := [ tup : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
        for tup in Fac do
            K := ExtendRelativeSplittingField(tup[1]);
        end for;
    end for;
else
    for f in fs do
        vprintf EndoFind : "Checking if %o has a root in %o\n", f, K;
        test := HasRoot(f, K);
        vprintf EndoFind : "done\n";
        if not test then
            Fac := [ tup : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
            for tup in Fac do
                K := ExtendRelativeSplittingField(tup[1]);
            end for;
        end if;
    end for;
end if;
return K;

end intrinsic;


intrinsic RelativeSplittingField(f::RngUPolElt) -> FldNum
{Returns a splitting field of the polynomials f over its base ring.}

return RelativeSplittingField([ f ]);

end intrinsic;


intrinsic RelativeSplittingFieldExtra(fs::SeqEnum) -> FldNum
{Returns a splitting field of the polynomials fs over their common base ring
together with an infinite place.}

K := RelativeSplittingField(fs);
SetInfinitePlaceUpwards(K);
return K;

end intrinsic;


intrinsic RelativeSplittingFieldExtra(f::RngUPolElt) -> FldNum
{Returns a splitting field of the polynomial f over its base ring together
with an infinite place.}

return RelativeSplittingFieldExtra([ f ]);

end intrinsic;


intrinsic RelativeCompositum(K::Fld, L::Fld) -> Fld
{Returns the compositum of the fields K and L over their common base ring, by
taking the splitting field of the polynomial defining L to extend K.}

if IsQQ(K) and IsQQ(L) then
    M := K;
    phiK := hom<K -> M | >;
    phiL := hom<L -> M | >;
elif IsQQ(K) then
    M := L;
    phiK := hom<K -> M | >;
    phiL := hom<L -> M | L.1>;
elif IsQQ(L) then
    M := K;
    phiK := hom<K -> M | K.1>;
    phiL := hom<L -> M | >;
else

    F := BaseRing(K); M := K;
    R<x> := PolynomialRing(F);
    g := MinimalPolynomial(L.1, F);
    tup := Factorization(g, K)[1];
    M := ExtendRelativeSplittingField(tup[1]);
    M := ImproveField(M);
    testK, phiK := IsSubfield(K, M); testL, phiL := IsSubfield(L, M);
end if;
return M, [* phiK, phiL *];

end intrinsic;


intrinsic RelativeCompositum(Ks::List) -> Fld, List
{Returns the compositum of the fields in Ks over their common base ring,
recursively adjoining the splitting field of the last factor.}

if #Ks eq 1 then
    return Ks[1], hom< Ks[1] -> Ks[1] | Ks[1].1 >;
elif #Ks eq 2 then
    L, psis := RelativeCompositum(Ks[1], Ks[2]);
    return L, psis;
end if;
L, psis := RelativeCompositum(Ks[1..(#Ks - 1)]);
M, phis := RelativeCompositum(L, Ks[#Ks]);
return M, [* psi * phis[1] : psi in psis *] cat [* phis[2] *];

end intrinsic;


//intrinsic RelativeFixedField(L::Fld, gens::SeqEnum) -> Fld
//{Returns the fixed subfield of L under the automorphisms in gens, considered as
//a field over the base ring of L.}
//
///* TODO: FixedField itself seems to work */
//dL := Degree(L); F := BaseRing(L);
//Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
//Ker := &meet[ Kernel(M) : M in Ms ];
///* Get F-basis of fixed field over F */
//B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
///* Subfield automatically creates K as extension of F */
//K := sub< L | B >;
//return K;
//
//end intrinsic;


intrinsic FixedFieldExtra(L::Fld, gens::SeqEnum) -> Fld
{Returns the fixed subfield of L under the automorphisms in gens, considered as
a field over the base ring of L.}

K := FixedField(L, gens);
SetInfinitePlaceDownwards(L, K);
return K;

end intrinsic;


intrinsic TransferInfinitePlace(h::Map, iotaK::.) -> .
{Return the image of the infinite place iotaK under the isomorphism h.}

/* TODO: This precision should not be global */
prec := 100;
K := Domain(h); L := Codomain(h);
if IsQQ(Domain(h)) then
    return InfinitePlaces(L)[1];
end if;
for iotaL in InfinitePlaces(L) do
    evL1 := Evaluate(h(K.1), iotaL : Precision := prec);
    evL2 := ComplexConjugate(evL1);
    evK := Evaluate(K.1, iotaK : Precision := prec);
    if Abs(evL1 - evK) lt 10^(-prec + 10) or Abs(evL2 - evK) lt 10^(-prec + 10) then
        return iotaL;
    end if;
end for;

end intrinsic;


intrinsic IsIsomorphicExtra(K::Fld, L::Fld) -> BoolElt, Map
{Returns the usual isomorphism and transfer the infinite place from K to L.}

if (IsQQ(K) and not IsQQ(L)) or (IsQQ(L) and not IsQQ(K)) then
    return false;
end if;
if IsQQ(K) and IsQQ(L) then
    return true, hom< Rationals() -> Rationals() | >;
end if;
test, h := IsIsomorphic(K, L);
if test then
    L`iota := TransferInfinitePlace(h, K`iota);
end if;
return test, h;

end intrinsic;


intrinsic TransferMatrices(As::SeqEnum, K::Fld, L::Fld) -> SeqEnum
{Returns the matrices in As transformed by some isomorphism from K to L.}

test, h := IsIsomorphicExtra(K, L);
return [ Matrix([ [ h(c) : c in Eltseq(row) ] : row in Rows(A) ]) : A in As ];

end intrinsic;
