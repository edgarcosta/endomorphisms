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


intrinsic FieldDescription(K::Fld) -> SeqEnum
{Returns a list describing the field K.}

/* Less structured fields encountered in endomorphism lattice */
if not assigned K`base then
    return [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
end if;
/* Now fields with extra structure */
if IsQQ(K) then
    F := K;
else
    F := sub< K | K`base_gen >;
end if;
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
K := AbsoluteField(NumberField(f));
K`base := F;
K`base_gen := K ! F.1;
K`CC := F`CC;
h := CanonicalInclusionMap(F, K);
K`iota := AscendInfinitePlace(F, K, h);
return K;

end intrinsic;


intrinsic EmbedAtInfinitePlace(f::RngUPolElt) -> RngUPolElt
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


intrinsic EmbedAtInfinitePlace(f::RngMPolElt) -> RngMPolElt
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


intrinsic EmbedAtInfinitePlace(fs::SeqEnum) -> SeqEnum
{Returns the list of polynomials fs considered as complex polynomials to
precision prec.}

return [ EmbedAtInfinitePlace(f) : f in fs ];

end intrinsic;


intrinsic DescendInfinitePlace(L::Fld, K::Fld, h::Map) -> Fld
{Descends infinite place from L to domain K of h.}

assert K eq Domain(h);
assert L eq Codomain(h);
assert Precision(K`CC) eq Precision(L`CC);
CC := K`CC; genK := K.1; genL := h(K.1);
for iotaK in InfinitePlacesExtra(K) do
    genKCC := CC ! EvaluateExtra(genK, iotaK);
    genLCC := CC ! EvaluateExtra(genL, L`iota);
    if Abs(genKCC - genLCC) lt CC`epscomp then
        return iotaK;
    end if;
end for;

end intrinsic;


intrinsic AscendInfinitePlace(K::Fld, L::Fld, h::Map) -> Fld
{Descends infinite place from L to domain K of h.}

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

end intrinsic;



intrinsic CanonicalInclusionMap(K::Fld, L::Fld) -> Map
{Gives the canonical inclusion map from K into L.}

assert IsSubfield(K, L);
if IsQQ(K) then
    return hom< K -> L | >;
else
    return hom< K -> L | L ! K.1 >;
end if;

end intrinsic;


intrinsic ImproveFieldExtra(K::Fld) -> Fld, Map
{Polredbestabs plus attribute transfer.}

K0, hKK0 := Polredbestabs(K);
TransferAttributes(K, K0, hKK0);
return K0, hKK0;

end intrinsic;


intrinsic FixedFieldExtra(L::Fld, gens::SeqEnum) -> Fld
{Returns the fixed subfield K of L under the automorphisms in gens, along with the inclusion of K in L.}

dL := Degree(L);
Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
Ker := &meet[ Kernel(M) : M in Ms ];
B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
K := sub< L | B >; hKL := CanonicalInclusionMap(K, L);
K`base := L`base;
K`base_gen := K ! L`base_gen;
K`CC := L`CC;
K`iota := DescendInfinitePlace(L, K, hKL);
if Type(K) eq FldRat then
    return K, hom< K -> L | >;
end if;
K0, hKK0 := ImproveFieldExtra(K); hKK0i := Inverse(hKK0);
h := hom< K0 -> L | hKL(hKK0i(K0.1)) >;
return K0, h;

end intrinsic;


intrinsic TransferAttributes(K::Fld, L::Fld, h::Map)
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

/* Initiate and make sure that the field of K is good for comparison purposes */
if K eq F then
    K := F; K`base := K; K`base_gen := K.1;
end if;
assert Precision(Parent(aCCs[1])) ge Precision(K`CC);

/* Iterative extension */
tupsa := [ ];
for aCC in aCCs do
    vprint EndoFind, 2 : "Before extension:";
    vprint EndoFind, 2 : K;
    vprint EndoFind : "";
    K, tupsa := ExtendNumberFieldExtra(K, tupsa, aCC);
    vprint EndoFind, 2 : "";
    vprint EndoFind, 2 : "After extension:";
    vprint EndoFind, 2 : K;
    vprint EndoFind, 2 : "";
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
        vprint EndoFind, 2 : iotaL;
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


intrinsic SplittingFieldExtra(aCCs::SeqEnum[FldComElt], F::Fld) -> Fld, SeqEnum
{Given complex number aCCs and a NumberFieldExtra F, finds the splitting field generated by the elements of aCCs and realizes said elements in that field.}

/* Initiate and make sure that the field of K is good for comparison purposes */
K := F; K`base := K; K`base_gen := K.1;
assert Precision(Parent(aCCs[1])) ge Precision(K`CC);

/* Iterative extension */
tupsa := [ ];
for aCC in aCCs do
    vprint EndoFind, 2 : "Before extension:";
    vprint EndoFind, 2 : K;
    vprint EndoFind : "";
    K, tupsa := ExtendSplittingFieldExtra(K, tupsa, aCC);
    vprint EndoFind, 2 : "";
    vprint EndoFind, 2 : "After extension:";
    vprint EndoFind, 2 : K;
    vprint EndoFind, 2 : "";
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
        vprint EndoFind, 2 : iotaL;
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
            print tupa;
            print h(a);
            print aCC1;
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
        print test;

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
for aCC in aCCs do
    K, tupsa := ExtendNumberFieldExtra(K, tupsa, aCC : minpolQQ := gQQ);
end for;
return K, tupsa;

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


intrinsic CoerceToSubfieldMatrix(M::., L::Fld, K::Fld, h::Map) -> .
{Realizes M as a matrix over K.}

return Matrix([ [ CoerceToSubfieldElement(a, L, K, h) : a in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;
