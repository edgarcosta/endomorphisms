/***
 *  Number field functionality and attributes
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


declare attributes FldNum : base, base_gen, CC, iota;
declare attributes FldRat : base, base_gen, CC, iota;


intrinsic IsQQ(K::Fld) -> BoolElt
{Returns whether or not the field K equals QQ.}

return Type(K) eq FldRat;

end intrinsic;


intrinsic BaseFieldExtra(K::Fld) -> Fld, .
{Returns distinguished subfield and map to it.}

F := K`base;
if IsQQ(F) then
    FinK := F;
    h := hom< F -> FinK | >;
else
    FinK := sub< K | K`base_gen >;
    h := hom< F -> FinK | K`base_gen >;
end if;
return F, FinK, h;

end intrinsic;


intrinsic FieldDescription(K::Fld) -> SeqEnum
{Returns a list describing the field K.}

/* Less structured fields encountered in endomorphism lattice */
if not assigned K`base then
    return [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
end if;
/* Now fields with extra structure */
_, F, _ := BaseFieldExtra(K);
if IsQQ(F) then
    return [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
    return [ [ Rationals() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1, F)) ];
end if;

end intrinsic;


intrinsic ElementDescription(r::.) -> .
{Returns a list describing the field element r.}

K := Parent(r);
if IsQQ(K`base) then
    return [ Rationals() ! c : c in Eltseq(r) ];
else
    return [ [ Rationals() ! c : c in Eltseq(d) ] : d in Eltseq(r) ];
end if;

end intrinsic;


intrinsic InfinitePlacesExtra(K::Fld) -> SeqEnum
{The infinite places, but now in an actually useful form.}

return [ K`CC ! tup[1] : tup in Roots(MinimalPolynomial(K.1), ComplexField(Precision(K`CC) + 20)) ];

end intrinsic;


intrinsic EvaluateExtra(r::., iota::.) -> .
{Evaluates infinite place.}

seq := Eltseq(r);
return &+[ seq[i]*iota^(i - 1) : i in [1..#seq] ];

end intrinsic;


intrinsic EvaluateMatrixExtra(M::., iota::.) -> .
{Evaluates infinite place.}

return Matrix([ [ EvaluateExtra(r, iota) : r in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic RationalsExtra(prec::RngIntElt) -> FldNum
{Returns the number field defined by f along with an infinite place.}

K := Rationals();
K`base := K;
K`base_gen := K ! 1;
K`CC := ComplexFieldExtra(prec);
K`iota := InfinitePlacesExtra(K)[1];
return K;

end intrinsic;


intrinsic NumberFieldExtra(f::RngUPolElt) -> FldNum
{Returns the number field defined by f along with an infinite place.}

F := BaseRing(f);
K<r> := NumberField(f);
K := AbsoluteField(K);
K`base := F;
K`base_gen := K ! F.1;
K`CC := F`CC;
h := CanonicalInclusionMap(F, K);
K`iota := AscendInfinitePlace(F, K, h);
K0, hKK0 := ImproveFieldExtra(K);
return K0, hKK0(K ! r), hKK0(K ! F.1);

end intrinsic;


intrinsic EmbedAtInfinitePlacePolynomial(f::RngUPolElt) -> RngUPolElt
{Returns the polynomial f considered as a complex polynomial to precision
prec.}

K := BaseRing(f);
RCC := PolynomialRing(K`CC);
if IsZero(f) then
    return RCC ! 0;
else
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ EvaluateExtra(MonomialCoefficient(f, mon), K`iota) * RCC.1^Degree(mon) : mon in mons ];
end if;

end intrinsic;


intrinsic EmbedAtInfinitePlacePolynomial(f::RngMPolElt) -> RngMPolElt
{Returns the polynomial f considered as a complex polynomial to precision
prec.}

K := BaseRing(Parent(f));
RCC := PolynomialRing(K`CC, #GeneratorsSequence(Parent(f)));
if IsZero(f) then
    return RCC ! 0;
else
    prec := Precision(BaseRing(RCC));
    mons := Monomials(f);
    return &+[ EvaluateExtra(MonomialCoefficient(f, mon), K`iota) * Monomial(RCC, Exponents(mon)) : mon in mons ];
end if;

end intrinsic;


intrinsic EmbedAtInfinitePlacePolynomials(fs::SeqEnum) -> SeqEnum
{Returns the list of polynomials fs considered as complex polynomials to
precision prec.}

return [ EmbedAtInfinitePlacePolynomial(f) : f in fs ];

end intrinsic;


intrinsic DescendInfinitePlace(L::Fld, K::Fld, h::Map) -> Fld
{Descends infinite place from codomain L to domain K of h.}

assert K eq Domain(h);
assert L eq Codomain(h);
assert Precision(K`CC) eq Precision(L`CC);
CC := K`CC; genK := K.1; genL := h(K.1);
for iotaK in InfinitePlacesExtra(K) do
    genKCC := CC ! EvaluateExtra(genK, iotaK);
    genLCC := CC ! EvaluateExtra(genL, L`iota);
    //vprint EndoFind, 3 : RealField(5) ! Abs(genKCC - genLCC);
    if Abs(genKCC - genLCC) lt CC`epscomp then
        return iotaK;
    end if;
end for;
error "Failed to descend infinite place";

end intrinsic;


intrinsic AscendInfinitePlace(K::Fld, L::Fld, h::Map) -> Fld
{Ascends infinite place from domain K to codomain L of h.}

assert K eq Domain(h);
assert L eq Codomain(h);
assert Precision(K`CC) eq Precision(L`CC);
CC := K`CC; genK := K.1; genL := h(K.1);
for iotaL in InfinitePlacesExtra(L) do
    genKCC := CC ! EvaluateExtra(genK, K`iota);
    genLCC := CC ! EvaluateExtra(genL, iotaL);
    if Abs(genKCC - genLCC) lt CC`epscomp then
        return iotaL;
    end if;
end for;
error "Failed to ascend infinite place";

end intrinsic;


intrinsic CanonicalInclusionMap(K::Fld, L::Fld) -> Map
{Gives the canonical inclusion map from K into L.}

if IsQQ(K) then
    return hom< K -> L | >;
else
    assert IsSubfield(K, L);
    return hom< K -> L | L ! K.1 >;
end if;

end intrinsic;


intrinsic SubfieldExtra(L::Fld, seq::.) -> .
{Gives subfield of L generated by seq.}

if IsQQ(L) then
    K := L;
    h := hom< K -> L | >;
else
    K, h := sub< L | Eltseq(seq) cat [ L`base_gen ] >;
end if;
RestrictAttributesExtra(L, K, h);
return K, h;

end intrinsic;


intrinsic RestrictAttributesExtra(L::Fld, K::Fld, h::Map)
{Restrict the attributes from L to K using the homomorphism h.}

assert K eq Domain(h);
assert L eq Codomain(h);
K`base := L`base;
K`base_gen := CoerceToSubfieldElement(L`base_gen, L, K, h);
K`CC := L`CC; CC := K`CC;
if IsQQ(K) then
    K`iota := InfinitePlacesExtra(K)[1]; return;
end if;
for iotaK in InfinitePlacesExtra(K) do
    evL := EvaluateExtra(h(K.1), L`iota);
    evK := EvaluateExtra(K.1, iotaK);
    if Abs(evL - evK) lt CC`epscomp then
        K`iota := iotaK; return;
    end if;
end for;

end intrinsic;




intrinsic ImproveFieldExtra(K::Fld) -> Fld, Map
{Polredbestabs plus attribute transfer.}

K0, hKK0 := Polredbestabs(K);
TransferAttributesExtra(K, K0, hKK0);
return K0, hKK0;

end intrinsic;


intrinsic FixedFieldExtra(L::Fld, gens::SeqEnum) -> Fld
{Returns the fixed subfield K of L under the automorphisms in gens, along with the inclusion of K in L.}

if #gens eq 0 then
    return L;
end if;
dL := Degree(L);
Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
Ker := &meet[ Kernel(M) : M in Ms ];
B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
K := sub< L | B >; hKL := CanonicalInclusionMap(K, L);
K`base := L`base;
K`base_gen := K ! L`base_gen;
K`CC := L`CC;
if Type(K) eq FldRat then
    K`iota := InfinitePlacesExtra(K)[1];
    return K, hom< K -> L | >;
end if;
K0, hKK0 := Polredbestabs(K);
hKK0i := Inverse(hKK0);
hK0L := hom< K0 -> L | hKL(hKK0i(K0.1)) >;
K0`base := K`base;
K0`base_gen := hKK0(K`base_gen);
K0`CC := K`CC;
K0`iota := DescendInfinitePlace(L, K0, hK0L);
/* TODO: This line should go, but that leads to problems in Sage */
//K0`iota := InfinitePlacesExtra(K0)[1];
return K0, hK0L;

end intrinsic;


intrinsic FixedGroupExtra(L::Fld, K::Fld, AutL::.) -> .
{Extension of FixedGroup.}

Gp, Gf, Gphi := Explode(AutL);
if IsQQ(L) then
    return Gp;
else
    return sub< Gp | [ g : g in Gp | (Gphi(g))(K.1) eq K.1 ] >;
end if;

end intrinsic;


intrinsic TransferAttributesExtra(K::Fld, L::Fld, h::Map)
{Transfer the attributes from K to L using the homomorphism h.}

assert K eq Domain(h);
assert L eq Codomain(h);
L`base := K`base;
L`base_gen := h(K`base_gen);
L`CC := K`CC; CC := K`CC;
if IsQQ(K) then
    L`iota := InfinitePlacesExtra(L)[1]; return;
end if;
for iotaL in InfinitePlacesExtra(L) do
    evL := EvaluateExtra(h(K.1), iotaL);
    evK := EvaluateExtra(K.1, K`iota);
    if Abs(evL - evK) lt CC`epscomp then
        L`iota := iotaL; return;
    end if;
end for;

end intrinsic;


intrinsic ConjugatePolynomial(h::Map, f::RngUPolElt) -> RngUPolElt
{Returns the transformation of the univariate polynomial f by the map h.}

K := Domain(h); L := Codomain(h); S := PolynomialRing(L);
return &+[ h(Coefficient(f, i))*S.1^i : i in [0..Degree(f)] ];

end intrinsic;


intrinsic ConjugateMatrix(h::Map, M::.) -> .
{Returns the transformation of the matrix M by the map h.}

return Matrix([ [ h(elt) : elt in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic NumberFieldExtra(aCCs::SeqEnum[FldComElt], F::Fld : K := F) -> Fld, SeqEnum
{Given complex number aCCs and a NumberFieldExtra F, finds the number field generated by the elements of aCCs and realizes said elements in that field.}

if #aCCs eq 0 then
    return F, [ ];
end if;

/* Initiate and make sure that the field of K is good for comparison purposes */
if K eq F then
    K := F; K`base := K; K`base_gen := K.1;
end if;
assert Precision(Parent(aCCs[1])) ge Precision(K`CC);

/* Iterative extension */
tupsa := [ ];
for aCC in aCCs do
    Knew, tupsa := ExtendNumberFieldExtra(K, tupsa, aCC);
    if Knew ne K then
        vprint EndoFind : "";
        vprint EndoFind : "After extension:";
        vprint EndoFind : Knew;
        vprint EndoFind : "";
    end if;
    K := Knew;
end for;

/* Sanity check before returning */
F := K`base; CC := K`CC;
genFCC0 := CC ! EvaluateExtra(F.1, F`iota);
genFCC1 := CC ! EvaluateExtra(K`base_gen, K`iota);
assert Abs(genFCC1 - genFCC0) lt CC`epscomp;
for tupa in tupsa do
    aCC0 := CC ! tupa[2];
    aCC1 := CC ! EvaluateExtra(tupa[1], K`iota);
    assert Abs(aCC1 - aCC0) lt CC`epscomp;
end for;
return K, [ tupa[1] : tupa in tupsa ];

end intrinsic;


intrinsic ExtendNumberFieldExtra(K::Fld, tupsa::SeqEnum, anewCC::FldComElt : minpolQQ := 0) -> Fld, SeqEnum
{Let K be a NumberFieldExtra with elements tupsa, and let anewCC be a complex number. This function returns the field generated by K and anewCC, and transports tupsa to that field.}

/* Determine minimal polynomial over K */
gK, gQQ := MinimalPolynomialExtra(anewCC, K : minpolQQ := minpolQQ);
if Degree(gK) eq 1 then
    anew := -Coefficient(gK, 0)/Coefficient(gK, 1);
    Append(~tupsa, <anew, anewCC>);
    return K, tupsa;
end if;

/* Information about the base needed later */
F := K`base; genFCC0 := EvaluateExtra(F.1, F`iota);
CC := K`CC;

/* Get absolute field and we need an iso that respects results so far */
Lrel := NumberField(gK);
Labs := AbsoluteField(Lrel);
L, h := Polredbestabs(Labs);
//L := Labs; h := hom< L -> L | L.1 >;
anew := h(Labs ! Lrel.1);
f := MinimalPolynomial(K.1);
rtsf := RootsPari(f, L);
L`CC := CC; L`base := F;

/* Choose compatible root */
for rtf in rtsf do
    if IsQQ(K) then
        h := hom< K -> L | >;
    else
        h := hom< K -> L | rtf >;
    end if;
    L`base_gen := h(K`base_gen);

    for iotaL in InfinitePlacesExtra(L) do
        test := true;

        /* First test: generator of F */
        genFCC1 := EvaluateExtra(L`base_gen, iotaL);
        if not Abs(genFCC1 - genFCC0) lt CC`epscomp then
            test := false;
        end if;

        /* Second test: old a */
        for tupa in tupsa do
            a := tupa[1]; aCC0 := tupa[2];
            aCC1 := EvaluateExtra(h(a), iotaL);
            if not Abs(aCC1 - CC ! aCC0) lt CC`epscomp then
                test := false;
                break;
            end if;
        end for;

        /* Third test: new a */
        anewCC1 := EvaluateExtra(anew, iotaL);
        if not Abs(anewCC1 - CC ! anewCC) lt CC`epscomp then
            test := false;
        end if;

        if test then
            L`iota := iotaL;
            newa := < anew, anewCC >;
            tupsa := [ < h(tupa[1]), tupa[2] > : tupa in tupsa ] cat [ newa ];
            return L, tupsa;
        end if;
    end for;
end for;
error "Failed to find a relative number field";

end intrinsic;


intrinsic ExtendNumberFieldExtra(g::RngUPolElt) -> .
{Let g be a polynomial over the NumberFieldExtra K. This function returns the
field generated by g with the same base, plus an inclusion from K into the new field.}

K := BaseRing(Parent(g));
Lrel := NumberField(g);
Labs := AbsoluteField(Lrel);
L, h := Polredbestabs(Labs);

f := MinimalPolynomial(K.1);
rtsf := RootsPari(f, L);
L`CC := K`CC; L`base := K`base;
if IsQQ(K) then
    h := hom< K -> L | >;
else
    h := hom< K -> L | rtsf[1] >;
end if;
L`base_gen := h(K`base_gen);
L`iota := AscendInfinitePlace(K, L, h);
return L, h;

end intrinsic;


intrinsic SplittingFieldExtra(aCCs::SeqEnum[FldComElt], F::Fld) -> Fld, SeqEnum
{Given complex number aCCs and a NumberFieldExtra F, finds the splitting field generated by the elements of aCCs and realizes said elements in that field.}

if #aCCs eq 0 then
    return F, [ ];
end if;

/* Initiate and make sure that the field of K is good for comparison purposes */
K := F; K`base := K; K`base_gen := K.1;
assert Precision(Parent(aCCs[1])) ge Precision(K`CC);

/* Iterative extension */
tupsa := [ ];
for aCC in aCCs do
    Knew, tupsa := ExtendSplittingFieldExtra(K, tupsa, aCC);
    if Knew ne K then
        vprint EndoFind : "";
        vprint EndoFind : "After extension:";
        vprint EndoFind : Knew;
        vprint EndoFind : "";
    end if;
    K := Knew;
end for;

/* Sanity check before returning */
F := K`base; CC := K`CC;
genFCC0 := CC ! EvaluateExtra(F.1, F`iota);
genFCC1 := CC ! EvaluateExtra(K`base_gen, K`iota);
assert Abs(genFCC1 - genFCC0) lt CC`epscomp;
for tupa in tupsa do
    aCC0 := CC ! tupa[2];
    aCC1 := CC ! EvaluateExtra(tupa[1], K`iota);
    assert Abs(aCC1 - aCC0) lt CC`epscomp;
end for;
return K, [ tupa[1] : tupa in tupsa ];

end intrinsic;


intrinsic ExtendSplittingFieldExtra(K::Fld, tupsa::SeqEnum, anewCC::FldComElt) -> Fld, SeqEnum
{Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a complex number. This function returns the splitting field generated by K and anewCC, and transports tupsa to that field.}

if IsQQ(K`base) then
    return ExtendSplittingFieldExtraQQ(K, tupsa, anewCC);
end if;
return ExtendSplittingFieldExtraGen(K, tupsa, anewCC);

end intrinsic;


intrinsic ExtendSplittingFieldExtraQQ(K::Fld, tupsa::SeqEnum, anewCC::FldComElt : minpolQQ := 0) -> Fld, SeqEnum
{Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a complex number. This function returns the splitting field generated by K and anewCC, and transports tupsa to that field.}

/* Determine minimal polynomial over K */
gK, gQQ := MinimalPolynomialExtra(anewCC, K : minpolQQ := minpolQQ);
if Degree(gK) eq 1 then
    anew := -Coefficient(gK, 0)/Coefficient(gK, 1);
    Append(~tupsa, <anew, anewCC>);
    return K, tupsa;
end if;

/* Information about the base needed later */
F := K`base; genFCC0 := EvaluateExtra(F.1, F`iota);
CC := K`CC;

/* Get absolute field and we need an iso that respects results so far */
Lrel := Compositum(K, SplittingFieldPari(gQQ));
Labs := AbsoluteField(Lrel);
L, h := Polredbestabs(Labs);
rtsg := RootsPari(gQQ, L);
f := MinimalPolynomial(K.1);
rtsf := RootsPari(f, L);
L`CC := CC; L`base := F;

/* Choose compatible root */
for rtf in rtsf do
    if IsQQ(K) then
        h := hom< K -> L | >;
    else
        h := hom< K -> L | rtf >;
    end if;
    L`base_gen := h(K`base_gen);

    for iotaL in InfinitePlacesExtra(L) do
        test := true;

        /* First test: generator of F */
        genFCC1 := EvaluateExtra(L`base_gen, iotaL);
        if not Abs(genFCC1 - genFCC0) lt CC`epscomp then
            test := false;
        end if;

        /* Second test: old a */
        for tupa in tupsa do
            a := tupa[1]; aCC0 := tupa[2];
            aCC1 := EvaluateExtra(h(a), iotaL);
            if not Abs(aCC1 - CC ! aCC0) lt CC`epscomp then
                test := false;
                break;
            end if;
        end for;

        /* Third test: new a */
        if test then
            test := false;
            for rtg in rtsg do
                anewCC1 := EvaluateExtra(rtg, iotaL);
                if Abs(anewCC1 - CC ! anewCC) lt CC`epscomp then
                    test := true;
                    anew := rtg;
                end if;
            end for;
        end if;

        if test then
            L`iota := iotaL;
            newa := < anew, anewCC >;
            tupsa := [ < h(tupa[1]), tupa[2] > : tupa in tupsa ] cat [ newa ];
            return L, tupsa;
        end if;
    end for;
end for;
error "Failed to find a relative number field";

end intrinsic;


intrinsic ExtendSplittingFieldExtraGen(K::Fld, tupsa::SeqEnum, anewCC::FldComElt) -> Fld, SeqEnum
{Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a complex number. This function returns the splitting field generated by K and anewCC, and transports tupsa to that field.}

gK, gQQ := MinimalPolynomialExtra(anewCC, K);
CC := Parent(anewCC); aCCs := [ tup[1] : tup in Roots(gQQ, CC) ];
tupsaadd := [ ]; tupsanew := tupsa;
for aCC in aCCs do
    K, tupsanew := ExtendNumberFieldExtra(K, tupsanew, aCC : minpolQQ := gQQ);
    Append(~tupsaadd, tupsanew[#tupsanew]);
end for;
for tupa in tupsaadd do
    if Abs(tupa[2] - anewCC) lt CC`epscomp then
        return K, tupsanew[1..#tupsa] cat [ tupa ];
    end if;
end for;
error "No suitable root found in ExtendSplittingFieldExtraGen";

end intrinsic;


intrinsic CoerceToSubfieldElement(a::., L::Fld, K::Fld, h::Map) -> .
{Realizes a as an element of K.}

assert K eq Domain(h);
assert L eq Codomain(h);
if IsQQ(K) then
    return K ! a;
else
    return a @@ h;
end if;

end intrinsic;


intrinsic CoerceToSubfieldPolynomial(f::., L::Fld, K::Fld, h::Map) -> .
{Realizes f as a polynomial over K.}

assert K eq Domain(h);
assert L eq Codomain(h);
if Type(f) eq RngUPolElt then
    coeffs := Coefficients(f);
    coeffs0 := [CoerceToSubfieldElement(c, L, K, h) : c in coeffs ];
    return Polynomial(coeffs0);
elif Type(f) eq RngMPolElt then
    mons := Monomials(f);
    R := PolynomialRing(K, #GeneratorsSequence(Parent(f)));
    if IsZero(f) then
        return R ! 0;
    else
        return &+[ CoerceToSubfieldElement(MonomialCoefficient(f, mon), L, K, h) * Monomial(R, Exponent(mon)) : mon in mons ];
    end if;
end if;

end intrinsic;


intrinsic CoerceToSubfieldMatrix(M::., L::Fld, K::Fld, h::Map) -> .
{Realizes M as a matrix over K.}

return Matrix([ [ CoerceToSubfieldElement(a, L, K, h) : a in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic CompositumExtra(K::Fld, L::Fld) -> Fld
{Returns the compositum of the fields K and L over their common base field, by
taking the splitting field of the polynomial defining L to extend K. We require
L to be normal over the base.}
/* TODO: Outsource factorization and roots to more general Pari functions */

if IsQQ(K) and IsQQ(L) then
    M := K;
    hKM := hom<K -> M | >;
    hLM := hom<L -> M | >;
    return M, [* hKM, hLM *];
elif IsQQ(K) then
    M := L;
    hKM := hom<K -> M | >;
    hLM := hom<L -> M | L.1>;
    return M, [* hKM, hLM *];
elif IsQQ(L) then
    M := K;
    hKM := hom<K -> M | K.1>;
    hLM := hom<L -> M | >;
    return M, [* hKM, hLM *];
end if;

F, FinK, hFK := BaseFieldExtra(K);
F, FinL, hFL := BaseFieldExtra(L);

/* Find minimal polynomial of both fields over the base */
fK := MinimalPolynomial(K.1, FinK);
fL := MinimalPolynomial(L.1, FinL);
fLF := CoerceToSubfieldPolynomial(fL, FinL, F, hFL);
fKF := CoerceToSubfieldPolynomial(fK, FinK, F, hFK);

/* Take corresponding root and remember old one as well */
fLK := ConjugatePolynomial(hFK, fLF);
g := Factorization(fLK, K)[1][1];
M, rL, rK := NumberFieldExtra(g);
hKM := hom< K -> M | rK >;
hLM := hom< L -> M | rL >;

/* Move attributes */
M`CC := K`CC; M`base := K`base;
M`base_gen := hKM(K`base_gen);
M`iota := AscendInfinitePlace(K, M, hKM);
return M, [* hKM, hLM *];

end intrinsic;


intrinsic CompositumExtra(Ks::List) -> Fld, List
{Returns the compositum of the fields in Ks over their common base ring.}

if #Ks eq 1 then
    return Ks[1], hom< Ks[1] -> Ks[1] | Ks[1].1 >;
elif #Ks eq 2 then
    return CompositumExtra(Ks[1], Ks[2]);
end if;
L, psis := CompositumExtra(Ks[1..(#Ks - 1)]);
M, phis := CompositumExtra(L, Ks[#Ks]);
return M, [* psi * phis[1] : psi in psis *] cat [* phis[2] *];

end intrinsic;
