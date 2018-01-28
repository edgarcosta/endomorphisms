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

/* TODO: Simplify as far as possible in the case of doubly layered fields */

intrinsic IsQQ(F::Fld) -> BoolElt
{Returns whether or not the field F equals QQ.}

return Type(F) eq FldRat;

end intrinsic;


intrinsic HasBaseQQ(F::Fld) -> BoolElt
{Returns whether or not the field F has QQ as a base ring.}

return Type(BaseRing(F)) eq FldRat;

end intrinsic;


intrinsic GiveName(~K::Fld, F::Fld, str::MonStgElt)
{Give the name str to a generator of K over F if such a generator is indeed
present.}

if IsRelativeExtension(K, F) then
    AssignNames(~K, [ str ]);
end if;

end intrinsic;


intrinsic IsRelativeExtension(K::Fld, F::Fld) -> Fld
{Returns whether K is an extension of F. Equalities are not considered as
extensions.}

/* Over QQ, the field QQ is not an extension */
testrat := IsQQ(F) and IsQQ(K);
/* Over general fields, F is not an extension, and neither is a field with a different base */
testnonrat := (not IsQQ(F)) and (BaseRing(K) ne F);
/* Everything else is */
return not (testrat or testnonrat);

end intrinsic;


intrinsic MakeRelative(K::Fld, F::Fld) -> Fld
{Returns K as an extension of F, taking a linear extension if K equals F.}

if not IsRelativeExtension(K, F) then
    R<x> := PolynomialRing(F);
    K := NumberField(x - 1: DoLinearExtension := true);
end if;
return K;

end intrinsic;


intrinsic FieldDescription(K::Fld, F::Fld) -> SeqEnum
{Returns a list describing the field K as an extension of F.}

/* Make relative to deal with the case of a linear extension */
K := MakeRelative(K, F);
/* Use fewer brackets when the base field is QQ */
if IsQQ(F) then
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
    K_seq := [ [ Integers() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1)) ];
end if;
return K_seq;

end intrinsic;


intrinsic ElementDescription(r::., F::Fld) -> .
{Returns a list describing the field element r.}

K := Parent(r);
/* Use fewer brackets when possible in the next two cases */
if IsQQ(K) then
    return r;
elif (IsQQ(F) and (Degree(K) gt 1)) or (not IsRelativeExtension(K, F)) then
    return [ Rationals() ! c : c in Eltseq(r) ];
else
    return [ [ Rationals() ! c : c in Eltseq(d) ] : d in Eltseq(r) ];
end if;

end intrinsic;


intrinsic SetInfinitePlace(K::Fld, iota::.)
{Sets the infinite place of K to equal iota. This change is passed on
compatibly to the tower of fields to which K belongs.}

K`iota := iota;
if Type(K) eq FldRat then
    return;
end if;
F := BaseRing(K);
/* The next clause is clumsy, but needed because of problems when encountering
 * the rational field when it is left out */
if Type(F) eq FldRat then
    F`iota := InfinitePlaces(F)[1];
    return;
end if;
for iotaF in InfinitePlaces(F) do
    if Extends(iota, iotaF) then
        SetInfinitePlace(F, iotaF);
    end if;
end for;

end intrinsic;


intrinsic DefineOrExtendInfinitePlace(K::Fld)
{Assigns an infinite place to K that is compatible with the infinite place of
the base ring of K.}

F := BaseRing(K);
if not assigned F`iota or IsQQ(F) then
    SetInfinitePlace(K, InfinitePlaces(K)[1]);
else
    for iotaK in InfinitePlaces(K) do
        if Extends(iotaK, F`iota) then
            K`iota := iotaK;
        end if;
    end for;
end if;

end intrinsic;


intrinsic DefineOrExtendInfinitePlaceFunction(K::Fld) -> Fld
{Assigns an infinite place to K that is compatible with the infinite place of
the base ring F of K, then returns K. Used when invoking Sage.}

DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic RestrictInfinitePlace(L::Fld, K::Fld)
{Assigns an infinite place to K that is compatible with the infinite place of
the extension L of K. Both L and K should have the same base ring.}

if IsQQ(K) then
    SetInfinitePlace(K, InfinitePlaces(K)[1]); return;
else
    /* TODO: This invocation is a bit strange, see if it could give rise to
     * bugs */
    IsSubfield(K, L);
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
DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic NumberFieldExtra(f::RngUPolElt) -> FldNum
{Returns the number field defined by f along with an infinite place that
extends the place of the base field, if it has any.}

K := NumberField(f);
DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic EmbedAtInfinitePlace(f::RngUPolElt, RCC::RngUPol) -> RngUPolElt
{Returns the polynomial f considered as a complex polynomial to precision
prec.}

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
{Returns the list of polynomials fs considered as complex polynomials to
precision prec.}

R := Parent(fs[1]); d := #GeneratorsSequence(R);
F := BaseRing(R);
/* TODO: The field gets an infinite place if it did not have any, a slight
 * containment breach where side effects are concerned. One can be picky about
 * this in the future. */
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
{Returns a field isomorphic to K over its base ring such that the defining
polynomial of the new field is integral.}

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


intrinsic RelativeNumberField(K::Fld, F::Fld, f::RngUPolElt) -> Fld
{Given an extension K of F and a polynomial f over K, returns an extension of F
that contains both K and a root of f.}

if Degree(f) eq 1 then
    return K;
end if;

if K eq F then
    K := NumberField(f);
else
    /* This step can cost a lot of time */
    vprintf EndoFind : "Extending %o by root of %o over %o...\n", K, f, F;
    K := RelativeField(F, NumberField(f));
    vprintf EndoFind : "done\n";
end if;
K := ClearFieldDenominator(K);
return K;

end intrinsic;


intrinsic RelativeNumberFieldExtra(K::Fld, F::Fld, f::RngUPolElt) -> Fld
{Given an extension K of F and a polynomial f over K, returns an extension of F
that contains both K and a root of f.}

K := RelativeNumberFieldExtra(K, F, f);
DefineOrExtendInfinitePlace(K);
return K;

end intrinsic;


intrinsic ExtendRelativeSplittingField(K::Fld, F::Fld, f::RngUPolElt) -> Fld
{Given an extension K of F and a polynomial f over F, returns an extension of F
that contains both K and the splitting field of f over F.}

vprintf EndoFind : "Extending %o by splitting %o over %o\n...", K, f, F;
if Type(F) eq FldRat then
    K := SplittingField(f);
    vprintf EndoFind : "done\n";
    return K;
end if;

while true do
    factors := [ tup[1] : tup in Factorization(f, K) | Degree(tup[1]) gt 1 ];
    if #factors eq 0 then
        vprintf EndoFind : "done\n";
        return K;
    end if;
    K := RelativeNumberField(K, F, factors[1]);
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
{Returns a splitting field of the polynomials in fs over their common base
ring.}

F := BaseRing(fs[1]); K := F;
fs := Reverse(Sort(fs, CompareFields));
for f in fs do
    /* Note that the K in the next condition is updated as we proceed */
    vprintf EndoFind : "Checking if %o has a root in %o\n", f, K;
    test := HasRoot(f, K);
    vprintf EndoFind : "done\n";
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
{Returns a splitting field of the polynomials f over its base ring.}

return RelativeSplittingField([ f ]);

end intrinsic;


intrinsic RelativeSplittingFieldExtra(fs::SeqEnum) -> FldNum
{Returns a splitting field of the polynomials fs over their common base ring
together with an infinite place.}

K := RelativeSplittingField(fs);
DefineOrExtendInfinitePlace(K);
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
    M := ExtendRelativeSplittingField(K, F, tup[1]);
    M := ClearFieldDenominator(M);
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


intrinsic RelativeFixedField(L::Fld, gens::SeqEnum) -> Fld
{Returns the fixed subfield of L under the automorphisms in gens, considered as
a field over the base ring of L.}

/* TODO: FixedField itself seems to work */
dL := Degree(L); F := BaseRing(L);
Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
Ker := &meet[ Kernel(M) : M in Ms ];
/* Get F-basis of fixed field over F */
B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
/* Subfield automatically creates K as extension of F */
K := sub< L | B >;
return K;

end intrinsic;


intrinsic GeneralFixedField(L::Fld, gens::SeqEnum) -> Fld
{Returns the fixed subfield of L under the automorphisms in gens, considered as
a field over the base ring of L.}

/* TODO: This seems to have been fixed (pun not intended, darn) */
return FixedField(L, gens);
if HasBaseQQ(L) then
    return FixedField(L, gens), Rationals();
else
    return RelativeFixedField(L, gens), BaseRing(L);
end if;

end intrinsic;
