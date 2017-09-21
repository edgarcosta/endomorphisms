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


intrinsic IsQQ(F::Fld) -> BoolElt
{Whether or not F is QQ.}

return Type(F) eq FldRat;

end intrinsic;


intrinsic HasBaseQQ(F::Fld) -> BoolElt
{Whether or not F has QQ as a base ring.}

return Type(BaseRing(F)) eq FldRat;

end intrinsic;


intrinsic GiveName(~K::Fld, F::Fld, str::MonStgElt)
{Give a name to a generator if such a generator is indeed present.}

if IsRelativeExtension(K, F) then
    AssignNames(~K, [ str ]);
end if;

end intrinsic;


intrinsic IsRelativeExtension(K::Fld, F::Fld) -> Fld
{Checks if K is given as an extension of F.}

testrat := IsQQ(F) and IsQQ(K);
testnonrat := (not IsQQ(F)) and (BaseRing(K) ne F);
return not (testrat or testnonrat);

end intrinsic;


intrinsic MakeExtension(K::Fld, F::Fld) -> Fld
{Our conventions on field extensions.}

// TODO: Both versions should work in the end, but the first is far more
// user-friendly
return K;
return MakeRelative(K, F);

end intrinsic;


intrinsic MakeRelative(K::Fld, F::Fld) -> Fld
{Takes a linear extension if needed.}

if not IsRelativeExtension(K, F) then
    R<x> := PolynomialRing(F);
    K := NumberField(x - 1: DoLinearExtension := true);
end if;
return K;

end intrinsic;


intrinsic FieldDescription(K::Fld, F::Fld) -> SeqEnum
{Gives a list describing the field K as an extension of F.}

K := MakeRelative(K, F);
if IsQQ(F) then
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
    K_seq := [ [ Integers() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1)) ];
end if;
return K_seq;

end intrinsic;


intrinsic ElementDescription(r::., F::Fld) -> .
{Gives a list describing the field element r.}

K := Parent(r);
if IsQQ(K) then
    return r;
elif (IsQQ(F) and (Degree(K) gt 1)) or (not IsRelativeExtension(K, F)) then
    return [ Rationals() ! c : c in Eltseq(r) ];
else
    return [ [ Rationals() ! c : c in Eltseq(d) ] : d in Eltseq(r) ];
end if;

end intrinsic;


intrinsic SetInfinitePlace(K::FldNum, iota::.)
{Creates a complex field with some extra needed parameters.}

K`iota := iota;
F := BaseRing(K);
if not IsQQ(F) then
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


intrinsic RestrictInfinitePlace(L::Fld, K::Fld)
{Restricts an infinite place over a relative field extension.}

if IsQQ(K) then
    SetInfinitePlace(K, InfinitePlaces(K)[1]); return;
else
    /* TODO: This invocation is ridiculous */
    IsSubfield(K, L);
    for iotaK in InfinitePlaces(K) do
        if Extends(L`iota, iotaK) then
            K`iota := iotaK; return;
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
    return MakeExtension(K, F);
end if;
r := K.1;
coeffs := Coefficients(MinimalPolynomial(r, Rationals()));
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
return MakeExtension(K, F);

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
        // FIXME: This step can cost a lot of time
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
// dichotomies. However, Magma then takes much more time for some reason.
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
return MakeExtension(K, F);

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
{Relative compositum.}

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
{Fixed subfield of L determined by the automorphisms in gens.}

dL := Degree(L); F := BaseRing(L);
Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
Ker := &meet[ Kernel(M) : M in Ms ];
B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
K := sub< L | B >;
return MakeExtension(K, F);

end intrinsic;


intrinsic GeneralFixedField(L::Fld, gens::SeqEnum) -> Fld
{Fixed subfield of L determined by the automorphisms in gens.}

if HasBaseQQ(L) then
    return MakeExtension(FixedField(L, gens), Rationals());
else
    return MakeExtension(RelativeFixedField(L, gens), BaseRing(L));
end if;

end intrinsic;
