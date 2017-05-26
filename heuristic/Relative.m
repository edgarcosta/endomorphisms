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


/* We keep track of an infinite place, that can be modified if the need ever
 * arises. */

declare attributes FldNum : iota;
declare attributes FldRat : iota;


intrinsic SetInfinitePlace(K::FldNum, iota::.)
{Creates a complex field with some extra needed parameters.}

K`iota := iota;
F := BaseRing(K);
if Type(F) eq FldNum then
    for iotaF in InfinitePlaces(F) do
        if Extends(iota, iotaF) then
            SetInfinitePlace(F, iotaF);
        end if;
    end for;
else
    SetInfinitePlace(F, InfinitePlaces(F)[1]);
end if;

end intrinsic;


intrinsic SetInfinitePlace(K::FldRat, iota::.)
{Creates a complex field with some extra needed parameters.}

K`iota := iota;

end intrinsic;


intrinsic DefineOrExtendInfinitePlace(K::Fld)
{Extends an infinite place over a relative field extension.}

F := BaseRing(K);
if not assigned F`iota or Type(F) eq FldRat then
    SetInfinitePlace(K, InfinitePlaces(K)[1]);
else
    for iotaK in InfinitePlaces(K) do
        if Extends(iotaK, F`iota) then
            K`iota := iotaK;
        end if;
    end for;
end if;

end intrinsic;


intrinsic RestrictInfinitePlace(L::Fld, K::Fld)
{Restricts an infinite place over a relative field extension.}

if Type(K) eq FldRat then
    SetInfinitePlace(K, InfinitePlaces(K)[1]);
else
    for iotaK in InfinitePlaces(K) do
        if Extends(L`iota, iotaK) then
            K`iota := iotaK;
            break;
        end if;
    end for;
end if;

end intrinsic;


intrinsic DefineOrExtendInfinitePlaceFunction(K::Fld) -> Fld
{Extends an infinite place over a relative field extension.}

DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic NumberFieldExtra(f::RngUPolElt) -> FldNum
{Creates a number field with an embedding.}

K := NumberField(f);
DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic EmbedAtInfinitePlace(f::RngUPolElt, RCC::RngUPol) -> RngUPolElt
{Embeds the polynomial f in R into RCC via the infinite place iota.}

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
{Embeds the polynomial f in R into RCC via the infinite place iota.}

if IsZero(f) then
    return RCC ! 0;
else
    F := BaseRing(Parent(f));
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ Evaluate(MonomialCoefficient(f, mon), F`iota : Precision := prec) * Monomial(RCC, Exponents(mon)) : mon in mons ];
end if;

end intrinsic;


intrinsic EmbedAsComplexPolynomials(fs::SeqEnum, prec::RngIntElt) -> SeqEnum
{Embeds a list of polynomials fs as complex polynomials to precision prec.}

R := Parent(fs[1]); d := #GeneratorsSequence(R);
F := BaseRing(R);
if not assigned F`iota then
    SetInfinitePlace(F, InfinitePlaces(F)[1]);
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


intrinsic ClearFieldDenominator(K::Fld) -> Fld
{Simplifies the defining polynomial of a field to an integral version.}

F := BaseRing(K); d := Degree(K);
if d eq 1 then
    return F;
end if;
r := K.1;
coeffs := Coefficients(MinimalPolynomial(r, Rationals()));
dens := Reverse([ Denominator(coeff) : coeff in coeffs ]);
primes := &join[ Set([ tup[1] : tup in Factorization(den) ]) : den in dens | den ne 0 ];
if #primes eq 0 then
    common_den := 1;
else
    common_den := &*[ p^Maximum([ Ceiling(Valuation(dens[k], p)/k) : k in [1..#dens] | dens[k] ne 0 ]) : p in primes ];
end if;
f := MinimalPolynomial(common_den*r);
K := NumberField(f);
return K;

end intrinsic;


intrinsic ExtendRelativeSplittingField(K::Fld, F::Fld, f::RngUPolElt : Optimize := false) -> FldNum
{Extension step for relative splitting fields.}

if Degree(F) eq 1 then
    K := SplittingField(f);
    if Optimize then
        K := OptimizedRepresentation(K);
    end if;
    return K;
end if;
while true do
    factors := [ tup[1] : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
    if #factors eq 0 then
        return K;
    end if;
    if K eq F then
        K := NumberField(factors[1]);
    else
        // TODO: This step can cost a lot of time
        K := RelativeField(F, NumberField(factors[1]));
    end if;
end while;

end intrinsic;


function CompareFields(K1, K2);
// Input:   Two subfields, fields, or polynomials.
// Output:  A comparison function: field with smaller degrees are smaller.

if Degree(K1) lt Degree(K2) then
    return -1;
elif Degree(K1) eq Degree(K2) then
    return 0;
else
    return 1;
end if;

end function;


intrinsic RelativeSplittingField(fs::SeqEnum) -> FldNum
{Determines a relative splitting field of the polynomials in fs.}

// FIXME: This code below is relatively bad; I would prefer to first define K
// as a linear extension and then to update it. This removes certain
// dichotomies. However, Magma is not happy with it.
F := BaseRing(fs[1]); K := F;
fs := Reverse(Sort(fs, CompareFields));
for f in fs do
    // Note that the K in the next condition is updated as we proceed
    if not HasRoot(f, K) then
        for tup in Factorization(f, K) do
            K := ExtendRelativeSplittingField(K, F, tup[1]);
            K := ClearFieldDenominator(K);
        end for;
    end if;
end for;
return K;

end intrinsic;


intrinsic RelativeSplittingField(f::RngUPolElt) -> FldNum
{Determines a relative splitting field of the polynomials f.}

return RelativeSplittingField([ f ]);

end intrinsic;


intrinsic RelativeSplittingFieldExtra(fs::SeqEnum) -> FldNum
{Creates a relative splitting field field with an embedding.}

K := RelativeSplittingField(fs);
DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic RelativeSplittingFieldExtra(f::RngUPolElt) -> FldNum
{Creates a relative splitting field field with an embedding.}

return RelativeSplittingFieldExtra([ f ]);

end intrinsic;


intrinsic RelativeCompositum(K::Fld, L::Fld) -> Fld
{Relative compositum.}

if Type(K) eq FldRat and Type(L) eq FldRat then
    M := K;
    phiK := hom<K -> M | >;
    phiL := hom<L -> M | >;
elif Type(K) eq FldRat then
    M := L;
    phiK := hom<K -> M | >;
    phiL := hom<L -> M | L.1>;
elif Type(L) eq FldRat then
    M := K;
    phiK := hom<K -> M | K.1>;
    phiL := hom<L -> M | >;
else
    F := BaseRing(K); M := K;
    R<x> := PolynomialRing(F);
    g := MinimalPolynomial(L.1, F);
    tup := Factorization(g, K)[1];
    M := ExtendRelativeSplittingField(K, F, tup[1]);
    M := ClearFieldDenominator(M);
    testK, phiK := IsSubfield(K, M); testL, phiL := IsSubfield(L, M);
end if;
return M, phiK, phiL;

end intrinsic;


intrinsic RelativeCompositumExtra(K::Fld, L::Fld) -> Fld
{Relative compositum with an embedding. The latter will in general not be compatible.}

M := RelativeCompositum(K, L);
DefineOrExtendInfinitePlace(M);
return M;

end intrinsic;
