/***
 *  Number field functionality and attributes
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

/* Main idea: All fields have QQ as a base in Magma for simplifying
* manipulations, and the inclusion of the base is remembered by specifying an
* image of its generator */

declare attributes FldQuad : base, base_gen, CC, iota, aut;
declare attributes FldNum : base, base_gen, CC, iota, aut;
declare attributes FldRat : base, base_gen, CC, iota, aut;

forward ExtendNumberFieldExtraStep;
forward ExtendSplittingFieldExtraStep;
forward ExtendSplittingFieldExtraStepQQ;
forward ExtendSplittingFieldExtraStepGen;


intrinsic IsQQ(K::Fld) -> BoolElt
{Returns whether or not the field K equals QQ (with or without extra attributes).}

return Type(K) eq FldRat;

end intrinsic;


intrinsic InclusionOfBaseExtra(K::Fld) -> Fld, .
{Returns distinguished subfield of K and the map from it. Third and fourth return value map to image.}

F := K`base;
if IsQQ(F) then
  FinK := F; h := hom<F -> FinK | >; hFK := hom< F -> K | >;
else
  FinK := sub< K | K`base_gen >; h := hom<F -> FinK | K`base_gen>; hFK := hom<F -> K | K`base_gen>;
end if;
return F, hFK, FinK, h;

end intrinsic;


intrinsic FieldDescriptionExtra(K::Fld) -> SeqEnum
{Returns a list describing the field K.}

// Less structured fields encountered in endomorphism lattice
if not assigned K`base then
  return [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
end if;
// Now fields with extra structure
_, _, F := InclusionOfBaseExtra(K);
if IsQQ(F) then
  return [ Rationals() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
else
  return [ [ Rationals() ! c : c in Eltseq(F ! coeff) ] : coeff in Eltseq(MinimalPolynomial(K.1, F)) ];
end if;

end intrinsic;


intrinsic ElementDescriptionExtra(r::.) -> .
{Returns a list describing the field element r.}

K := Parent(r);
if IsQQ(K) then
  return r;
elif IsQQ(K`base) then
  return [ Rationals() ! c : c in Eltseq(r) ];
else
  return [ [ Rationals() ! c : c in Eltseq(d) ] : d in Eltseq(r) ];
end if;

end intrinsic;


intrinsic InfinitePlacesExtra(K::Fld) -> SeqEnum
{The infinite places of K, represented by the roots of the generator in the associated complex field. No identification of complex conjugate places takes place.}

return [ tup[1] : tup in Roots(MinimalPolynomial(K.1), ComplexFieldExtra(Precision(K`CC) + 10)) ];

end intrinsic;


intrinsic EmbedExtra(r::. : iota := 0) -> .
{Embeds the element r via the infinite place iota.}

F := Parent(r);
if Type(iota) eq RngIntElt then iota := F`iota; end if;
seq := Eltseq(Numerator(r));
return &+[ seq[i]*iota^(i - 1) : i in [1..#seq] ]/Denominator(r);

end intrinsic;


intrinsic EmbedMatrixExtra(M::. : iota := 0) -> .
{Embeds the matrix M via the infinite place iota.}

return Matrix([ [ EmbedExtra(r : iota := iota) : r in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic EmbedPolynomialExtra(f::RngUPolElt : iota := 0) -> RngUPolElt
{Embeds the polynomial f via the infinite place of its base.}

K := BaseRing(f);
RCC := PolynomialRing(Parent(K`iota));
if IsZero(f) then return RCC ! 0; end if;
mons := Monomials(f);
return &+[ EmbedExtra(MonomialCoefficient(f, mon) : iota := iota) * RCC.1^Degree(mon) : mon in mons ];

end intrinsic;


intrinsic EmbedPolynomialExtra(f::RngMPolElt) -> RngMPolElt
{Embeds the polynomial f via the infinite place of its base.}

K := BaseRing(Parent(f));
RCC := PolynomialRing(Parent(K`iota), #GeneratorsSequence(Parent(f)));
if IsZero(f) then return RCC ! 0; end if;
mons := Monomials(f);
return &+[ EmbedExtra(MonomialCoefficient(f, mon)) * Monomial(RCC, Exponents(mon)) : mon in mons ];

end intrinsic;


intrinsic EmbedPolynomialsExtra(fs::SeqEnum) -> SeqEnum
{Embeds the polynomials fs via the infinite place of their common base.}

return [ EmbedPolynomialExtra(f) : f in fs ];

end intrinsic;


intrinsic RationalsExtra(prec::RngIntElt) -> FldNum
{Returns the rationals with itself as base and an infinite place with the given precision.}
  K := Rationals();
  K`base := K;
  K`base_gen := K ! 1;
  K`CC := ComplexFieldExtra(prec);
  K`iota := InfinitePlacesExtra(K)[1];
return K;

end intrinsic;


intrinsic RationalsExtra() -> FldNum
{Default RationalsExtra with precision 100.}
  return RationalsExtra(100);
end intrinsic;

intrinsic RationalsExtra(b::BoolElt) -> FldNum
{ " } //"
  return RationalsExtra();
end intrinsic;

intrinsic BaseNumberFieldExtra(f::RngUPolElt, prec::RngIntElt : Simplify:=false) -> FldNum
{
  Returns the number field defined by f with itself as base and an infinite place with the given precision.
  The univariate polynomial f should be defined over QQ.
}

K := BaseRing(f);
assert IsQQ(K);
if Degree(f) eq 1 then
  L := RationalsExtra(prec);
  hKL := hom< K -> L | >;
  r := Roots(f, L)[1][1];
  return L, r, hKL;
end if;

Lrel<r> := NumberField(f);
L := AbsoluteField(Lrel);
L`base := L;
L`base_gen := L.1;
L`CC := ComplexFieldExtra(prec);
L`iota := InfinitePlacesExtra(L)[1];
hKL := hom< K -> L | >;

if Simplify then
  // Final improvement step before returning root
  L0, hLL0 := ImproveFieldExtra(L);
  return L0, hLL0(L ! r), hKL * hLL0;
end if;
return L, L!r, hKL;

end intrinsic;


intrinsic BaseNumberFieldExtra(f::RngUPolElt) -> FldNum
{Default BaseNumberFieldExtra defined by f with default precision.}

return BaseNumberFieldExtra(f, false);

end intrinsic;


intrinsic DescendInfinitePlace(L::Fld, K::Fld, h::Map) -> Fld
{Descends infinite place from the codomain L of h to its domain K.}

assert K eq Domain(h); assert L eq Codomain(h);
assert Precision(K`CC) eq Precision(L`CC);
CC := K`CC; genK := K.1; genL := h(K.1);
for iotaK in InfinitePlacesExtra(K) do
  genKCC := EmbedExtra(genK : iota := iotaK);
  genLCC := EmbedExtra(genL);
  if Abs(genKCC - genLCC) lt CC`epscomp then
    return iotaK;
  end if;
end for;
error "Failed to descend infinite place";

end intrinsic;


intrinsic AscendInfinitePlace(K::Fld, L::Fld, h::Map) -> Fld
{Ascends infinite place from the domain K of h to its codomain L.}

assert K eq Domain(h); assert L eq Codomain(h);
assert Precision(K`CC) eq Precision(L`CC);
CC := K`CC; genK := K.1; genL := h(K.1);

for iotaL in InfinitePlacesExtra(L) do
  genKCC := EmbedExtra(genK);
  genLCC := EmbedExtra(genL : iota := iotaL);
  if Abs(genKCC - genLCC) lt CC`epscomp then return iotaL; end if;
end for;
error "Failed to ascend infinite place";

end intrinsic;


intrinsic DescendAttributesExtra(L::Fld, K::Fld, h::Map)
{Descend the attributes from L to K using the homomorphism h.}

assert K eq Domain(h);
assert L eq Codomain(h);
K`base := L`base; K`base_gen := CoerceToSubfieldElement(L`base_gen, L, K, h);
K`CC := L`CC;
K`iota := DescendInfinitePlace(L, K, h);

end intrinsic;


intrinsic AscendAttributesExtra(L::Fld, K::Fld, h::Map)
{Ascend the attributes from K to L using the homomorphism h.}

assert K eq Domain(h);
assert L eq Codomain(h);
L`base := K`base; L`base_gen := h(K`base_gen);
L`CC := K`CC;
L`iota := AscendInfinitePlace(L, K, h);

end intrinsic;


intrinsic TransferAttributesExtra(K::Fld, L::Fld, h::Map)
{Transfer the attributes from K to L using the isomorphism h.}

assert K eq Domain(h);
assert L eq Codomain(h);
L`base := K`base;
L`base_gen := h(K`base_gen);
L`CC := K`CC; CC := K`CC;
if IsQQ(K) then L`iota := InfinitePlacesExtra(L)[1]; return; end if;

evK := EmbedExtra(K.1);
for iotaL in InfinitePlacesExtra(L) do
  evL := EmbedExtra(h(K.1) : iota := iotaL);
  if Abs(evL - evK) lt CC`epscomp then L`iota := iotaL; return; end if;
end for;

end intrinsic;


intrinsic CanonicalInclusionMap(K::Fld, L::Fld) -> Map
{Gives the canonical inclusion map from K into L, usually obtained by defining L as a NumberField over K.}

if IsQQ(K) then return hom< K -> L | >; end if;
assert IsSubfield(K, L);
return hom< K -> L | L ! K.1 >;

end intrinsic;


intrinsic SubfieldExtra(L::Fld, seq::. : Simplify:=true) -> .
{Gives the subfield of L generated by seq over its base, as a field plus an embedding plus the originals of the given sequence.}

if IsQQ(L) then
  K := L;
  h := hom< K -> L | >;
  return K, h, seq;
end if;

K := sub< L | seq cat [ L`base_gen ] >;
hKL := CanonicalInclusionMap(K, L);
if IsQQ(K) then
  DescendAttributesExtra(L, K, hKL);
  return K, hKL, [ K ! c : c in seq ];
end if;
if K eq L then
  return L, CanonicalInclusionMap(L, L), seq;
end if;

// Polishing
if Simplify then
  K0, hKK0 := Polredabs(K);
else
  K0 := K;
  hKK0 := IsQQ(K) select hom< K -> K0 | > else hom< K -> K0 | K0.1 >;
end if;
hKK0i := Inverse(hKK0);
hK0L := hom< K0 -> L | hKL(hKK0i(K0.1)) >;
DescendAttributesExtra(L, K0, hK0L);
return K0, hK0L, [ hKK0(K ! c) : c in seq ];

end intrinsic;


intrinsic ImproveFieldExtra(K::Fld) -> Fld, Map
{Polredabs plus attribute transfer. Returns the isomorphism.}

K0, hKK0 := Polredabs(K);
TransferAttributesExtra(K, K0, hKK0);
return K0, hKK0;

end intrinsic;


intrinsic FixedFieldExtra(L::Fld, gens::SeqEnum : Simplify:=true) -> Fld
{Returns the fixed subfield K of L over its base field under the automorphisms in gens. The inclusion of K in L is returned as a second value.}

if #gens eq 0 then return L, CanonicalInclusionMap(L, L); end if;
dL := Degree(L);
Ms := [ Matrix([ Eltseq(gen(L.1^i) - L.1^i) : i in [0..(dL - 1)] ]) : gen in gens ];
Ker := &meet[ Kernel(M) : M in Ms ];
B := [ &+[ b[i + 1]*L.1^i : i in [0..(dL - 1)] ] : b in Basis(Ker) ];
return SubfieldExtra(L, B : Simplify:=Simplify);

end intrinsic;


intrinsic FixedGroupExtra(L::Fld, K::Fld, h::.) -> .
{Returns the subgroup of the automorphism group of L over its base that fixes K elementwise.}

Gp, Gf, Gphi := AutomorphismGroupPari(L);
if IsQQ(L) then return Gp; end if;
r := h(K.1);
return sub< Gp | [ g : g in Gp | (Gphi(g))(r) eq r ] >;

end intrinsic;


intrinsic ConjugatePolynomial(h::Map, f::RngUPolElt) -> RngUPolElt
{Returns the transformation of the polynomial f by the map h.}

K := Domain(h); L := Codomain(h); S := PolynomialRing(L);
return &+[ h(Coefficient(f, i))*S.1^i : i in [0..Degree(f)] ];

end intrinsic;


intrinsic ConjugateMatrix(h::Map, M::.) -> .
{Returns the transformation of the matrix M by the map h.}

return Matrix([ [ h(elt) : elt in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic NumberFieldExtra(aCCs::SeqEnum[FldComElt], K::Fld : UpperBound:=16, DegreeDivides:=Infinity(), Simplify:=true) -> Fld, SeqEnum
{Given complex number aCCs and a number field K, finds the number field L over the base of K generated by the elements of aCCs and realizes said elements L. An inclusion map from K to L is also returned.}

// Trivial case
if #aCCs eq 0 then return K, [ ], CanonicalInclusionMap(K, K); end if;

// Initiate and make sure that the complex field of K is good for comparison purposes
L := K;
h := CanonicalInclusionMap(K, K);
assert Precision(Parent(aCCs[1])) ge Precision(L`CC);

// Iterative extension
tupsa := [ ];
for aCC in aCCs do
Lnew, tupsa, hnew := ExtendNumberFieldExtraStep(L, tupsa, aCC : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides, Simplify:=Simplify);
  if Lnew ne L then
    h := h*hnew;
    vprint EndoFind, 2 : "";
    vprint EndoFind, 2 : "Number field extended. Current field:";
    vprint EndoFind, 2 : Lnew;
  end if;
  L := Lnew;
end for;

// Sanity check before returning
F := L`base; CC := L`CC;
genKCC0 := EmbedExtra(K.1); genKCC1 := EmbedExtra(h(K.1));
assert Abs(genKCC1 - genKCC0) lt CC`epscomp;
for tupa in tupsa do
  aCC0 := tupa[2]; aCC1 := EmbedExtra(tupa[1]);
  assert Abs(aCC1 - aCC0) lt CC`epscomp;
end for;
return L, [ tupa[1] : tupa in tupsa ], h;

end intrinsic;


function ExtendNumberFieldExtraStep(K, tupsa, anewCC : UpperBound:=16, DegreeDivides:=Infinity(), Simplify:=true)
// Let K be a NumberFieldExtra with elements tupsa, and let anewCC be a complex
// number. This function returns the field generated by K and anewCC, and
// transports tupsa to that field. Keeps track of morphisms.

// Determine minimal polynomial over K
gK := MinimalPolynomialExtra(anewCC, K : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides);
// Done in case of degree 1
if Degree(gK) eq 1 then
  anew := -Coefficient(gK, 0)/Coefficient(gK, 1);
  Append(~tupsa, < anew, anewCC >);
  return K, tupsa, CanonicalInclusionMap(K, K);
end if;

// Information about the base
F := K`base; genFCC0 := EmbedExtra(F.1);
CC := K`CC; genKCC0 := EmbedExtra(K.1);

// Get absolute field and we need an iso that respects results so far
Lrel := NumberField(gK);
Labs := AbsoluteField(Lrel);
if Simplify then
  L, h := Polredabs(Labs);
else
  L := Labs;
  h := hom<Labs -> L | L.1>;
end if;
anew := h(Labs ! Lrel.1);
f := MinimalPolynomial(K.1);
rtsf := RootsPari(f, L);
L`CC := CC; L`base := F;

// Choose compatible root and embedding
for rtf in rtsf do
  h := IsQQ(K) select hom< K -> L | > else hom< K -> L | rtf >;
  L`base_gen := h(K`base_gen);
  for iotaL in InfinitePlacesExtra(L) do
    test := true;

    // First test: generator of K
    genKCC1 := EmbedExtra(h(K.1) : iota := iotaL);
    if not Abs(genKCC1 - genKCC0) lt CC`epscomp then
      test := false;
    end if;

    // Second test: old a
    for tupa in tupsa do
      a := tupa[1]; aCC0 := tupa[2];
      aCC1 := EmbedExtra(h(a) : iota := iotaL);
      if not Abs(aCC1 - aCC0) lt CC`epscomp then
        test := false;
        break;
      end if;
    end for;

    // Third test: new a
    anewCC1 := EmbedExtra(anew : iota := iotaL);
    if not Abs(anewCC1 - anewCC) lt CC`epscomp then
      test := false;
    end if;

    if test then
      L`iota := iotaL;
      newa := < anew, anewCC >;
      tupsa := [ < h(tupa[1]), tupa[2] > : tupa in tupsa ] cat [ newa ];
      return L, tupsa, h;
    end if;
    end for;
end for;
error "Failed to extend relative number field";

end function;


intrinsic NumberFieldExtra(f::RngUPolElt : prec:=false, Simplify:=false) -> .
{Given polynomial f, finds extension as NumberFieldExtra.}

K := BaseRing(f);
// if f in Z[x], cast it to Q[x]
if Type(K) eq RngInt then
  f := ChangeRing(f, Rationals());
  K := BaseRing(f);
end if;

// checking if we need call RationalsExtra(prec), and make f to have that base ring
if not assigned K`base or not assigned K`base`CC or prec cmpne false then
  // We deliberately ignore furnishing relative extensions... for now
  assert IsQQ(K);
  K := RationalsExtra(prec);
  R := PolynomialRing(K); f := R ! f;
  // since K already has precision, we can do false
  return NumberFieldExtra(f : prec:=false, Simplify:=Simplify);
end if;

if Degree(f) eq 1 then
  L := K;
  r := Roots(f)[1][1];
  hKL := CanonicalInclusionMap(K, L);
  return L, r, hKL;
end if;

Lrel<r> := NumberField(f);
L := AbsoluteField(Lrel);
L`base := K`base;
L`base_gen := L ! Lrel ! K`base_gen;
L`CC := K`CC;

// Inclusion map
hKL := IsQQ(K) select hom< K -> L | > else hom< K -> L | L ! Lrel ! K.1 >;
L`iota := AscendInfinitePlace(K, L, hKL);

if Simplify then
  // Final improvement step before returning
  L0, hLL0 := ImproveFieldExtra(L);
  return L0, hLL0(L ! r), hKL * hLL0;
end if;
return L, L!r, hKL;

end intrinsic;


intrinsic SplittingFieldExtra(f::RngUPolElt : prec:=false, Simplify:=false) -> .
{Given polynomial f, finds splitting field as NumberFieldExtra.}

K := BaseRing(f);
assert IsQQ(K);
// checking if we need call RationalsExtra(prec), and make f to have that base ring
if not assigned K`base or not assigned K`base`CC or prec cmpne false then
  // We deliberately ignore furnishing relative extensions... for now
  K := RationalsExtra(prec);
  R := PolynomialRing(K); f := R ! f;
  // since K already has precision, we can do false
  return SplittingFieldExtra(f : prec:=false, Simplify:=Simplify);
end if;

if Degree(f) eq 1 then
  L := K; r := Roots(f)[1][1]; hKL := CanonicalInclusionMap(K, L);
  return L, r, hKL;
end if;

Lrel<r>, rts := SplittingField(f);
L := AbsoluteField(Lrel);
L`base := K`base; L`base_gen := L ! Lrel ! K`base_gen; L`CC := K`CC;
rts := [ L ! rt : rt in rts ];

// Inclusion map
hKL := IsQQ(K) select hom< K -> L | > else hom< K -> L | L ! Lrel ! K.1 >;
L`iota := AscendInfinitePlace(K, L, hKL);

if Simplify then
  // Final improvement step before returning
  L0, hLL0 := ImproveFieldExtra(L);
  rts := [ hLL0(rt) : rt in rts ];
  return L0, hLL0(L ! r), hKL * hLL0, rts;
end if;
return L, L!r, hKL, rts;

end intrinsic;



intrinsic SplittingFieldExtra(aCCs::SeqEnum[FldComElt], K::Fld : UpperBound:=16, DegreeDivides:=Infinity(), Simplify:=true) -> Fld, SeqEnum
{Given complex number aCCs and a NumberFieldExtra K, finds the splitting field generated by the elements of aCCs and realizes said elements in that field. An inclusion map is also returned.}

// Trivial case
if #aCCs eq 0 then return K, [ ], CanonicalInclusionMap(K, K); end if;

// Initiate and make sure that the complex field of K is good for comparison purposes
L := K; h := CanonicalInclusionMap(K, K);
assert Precision(Parent(aCCs[1])) ge Precision(L`CC);

// Iterative extension
tupsa := [ ];
for aCC in aCCs do
  Lnew, tupsa, hnew := ExtendSplittingFieldExtraStep(L, tupsa, aCC : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides, Simplify:=Simplify);
  if Lnew ne L then
    h := h*hnew;
    vprint EndoFind, 2 : "";
    vprint EndoFind, 2 : "Number field extended. Current field:";
    vprint EndoFind, 2 : Lnew;
  end if;
  L := Lnew;
end for;

// Sanity check before returning
F := L`base;
CC := L`CC;
genKCC0 := EmbedExtra(K.1);
genKCC1 := EmbedExtra(h(K.1));
assert Abs(genKCC1 - genKCC0) lt CC`epscomp;
for tupa in tupsa do
  aCC0 := tupa[2];
  aCC1 := EmbedExtra(tupa[1]);
  assert Abs(aCC1 - aCC0) lt CC`epscomp;
end for;
return L, [ tupa[1] : tupa in tupsa ], h;

end intrinsic;


function ExtendSplittingFieldExtraStep(K, tupsa, anewCC : UpperBound:= 16, DegreeDivides:=Infinity(), Simplify:=true)
// Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a
// complex number. This function returns the splitting field generated by K and
// anewCC, and transports tupsa to that field. Keeps track of morphisms.

fn := IsQQ(K`base) select ExtendSplittingFieldExtraStepQQ else ExtendSplittingFieldExtraStepGen;
return fn(K, tupsa, anewCC : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides, Simplify:=Simplify);

end function;


function ExtendSplittingFieldExtraStepQQ(K, tupsa, anewCC : UpperBound:=16, DegreeDivides:=Infinity(), Simplify:=true)
// Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a
// complex number. This function returns the splitting field generated by K and
// anewCC, and transports tupsa to that field. Keeps track of morphisms.

// Determine minimal polynomial over K
//gK := MinimalPolynomialExtra(anewCC, K);
gK, gQQ := MinimalPolynomialExtra(anewCC, K : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides, UseQQ:=true);
// Done in case of degree 1
if Degree(gK) eq 1 then
  anew := -Coefficient(gK, 0)/Coefficient(gK, 1);
  Append(~tupsa, < anew, anewCC >);
  return K, tupsa, CanonicalInclusionMap(K, K);
end if;
//gK, gQQ := MinimalPolynomialExtra(anewCC, K : UseQQ := true);

// Information about the base needed later
F := K`base; genFCC0 := EmbedExtra(F.1); CC := K`CC;
genKCC0 := EmbedExtra(K.1);

// Get absolute field and we need an iso that respects results so far
Lrel := Compositum(K, SplittingFieldPari(gQQ));
Labs := AbsoluteField(Lrel);
if Simplify then
  L, h := Polredabs(Labs);
else
  L := Labs;
  h := hom<Labs -> L | L.1>;
end if;
L`CC := CC;
L`base := F;
rtsg := RootsPari(gQQ, L); f := MinimalPolynomial(K.1); rtsf := RootsPari(f, L);

// Choose compatible root
for rtf in rtsf do
  h := IsQQ(K) select hom< K -> L | > else hom< K -> L | rtf >;
  L`base_gen := h(K`base_gen);

  for iotaL in InfinitePlacesExtra(L) do
    test := true;

    // First test: generator of K
    genKCC1 := EmbedExtra(h(K.1) : iota := iotaL);
    if not Abs(genKCC1 - genKCC0) lt CC`epscomp then
      test := false;
    end if;

    // Second test: old a
    for tupa in tupsa do
      a := tupa[1]; aCC0 := tupa[2];
      aCC1 := EmbedExtra(h(a) : iota := iotaL);
      if not Abs(aCC1 - aCC0) lt CC`epscomp then
        test := false;
          break;
      end if;
    end for;

    // Third test: new a
    if test then
      test := false;
      for rtg in rtsg do
        anewCC1 := EmbedExtra(rtg : iota := iotaL);
        if Abs(anewCC1 - anewCC) lt CC`epscomp then
          test := true;
          anew := rtg;
        end if;
      end for;
    end if;

    if test then
      L`iota := iotaL;
      newa := < anew, anewCC >;
      tupsa := [ < h(tupa[1]), tupa[2] > : tupa in tupsa ] cat [ newa ];
      return L, tupsa, h;
    end if;
  end for;
end for;
error "Failed to extend relative splitting field";

end function;


function ExtendSplittingFieldExtraStepGen(K, tupsa, anewCC : UpperBound:=16, DegreeDivides:=Infinity(), Simplify:=true)
// Let K be a SplittingFieldExtra with elements tupsa, and let anewCC be a
// complex number. This function returns the splitting field generated by K and
// anewCC, and transports tupsa to that field. Keeps track of morphisms.

// Determine minimal polynomial over K
gK, gQQ := MinimalPolynomialExtra(anewCC, K : UpperBound:=UpperBound, DegreeDivides:=DegreeDivides);
// Done in case of degree 1
if Degree(gK) eq 1 then
  anew := -Coefficient(gK, 0)/Coefficient(gK, 1);
  Append(~tupsa, < anew, anewCC >);
  return K, tupsa, CanonicalInclusionMap(K, K);
end if;

// Information about the base needed later
F := K`base; CC := K`CC;
genFCC0 := EmbedExtra(F.1); genKCC0 := EmbedExtra(K.1);

// Get absolute field and we need an iso that respects results so far.
// We could also extend repeatedly, but this is more effective in practice.
// A more general Pari/GP version should be used
L, rtsg := SplittingField(gK);
h1 := IsQQ(K) select hom< K -> L | > else hom< K -> L | L ! (K.1) >;
if Simplify then
  L, h2 := Polredabs(L);
else
  h2 := IsQQ(L) select hom< L -> L | > else hom< L -> L | L.1 >;
end if;
h := h1*h2;
L`CC := CC;
L`base := F;
L`base_gen := h(K`base_gen);
anew := h2(rtsg[1]);

// Choose compatible embedding
for iotaL in InfinitePlacesExtra(L) do
  test := true;

  // First test: generator of K
  genKCC1 := EmbedExtra(h(K.1) : iota := iotaL);
  if not Abs(genKCC1 - genKCC0) lt CC`epscomp then
    test := false;
  end if;

  // Second test: old a
  for tupa in tupsa do
    a := tupa[1]; aCC0 := tupa[2];
    aCC1 := EmbedExtra(h(a) : iota := iotaL);
    if not Abs(aCC1 - aCC0) lt CC`epscomp then
      test := false;
      break;
    end if;
  end for;

  // Third test: new a
  if test then
    anewCC1 := EmbedExtra(anew : iota := iotaL);
    if not Abs(anewCC1 - anewCC) lt CC`epscomp then
      test := false;
    end if;
  end if;

  if test then
    L`iota := iotaL;
    newa := < anew, anewCC >;
    tupsa := [ < h(tupa[1]), tupa[2] > : tupa in tupsa ] cat [ newa ];
    return L, tupsa, h;
  end if;
end for;
error "Failed to extend relative splitting field";

end function;


intrinsic CoerceToSubfieldElement(a::., L::Fld, K::Fld, h::Map) -> .
{Coerces the element a to the subfield K of L via the map h from K to L.}

assert K eq Domain(h); assert L eq Codomain(h);
return IsQQ(K) select K ! a else K ! (a @@ h);

end intrinsic;


intrinsic CoerceToSubfieldPolynomial(f::., L::Fld, K::Fld, h::Map) -> .
{Coerces the polynomial f to the subfield K of L via the map h from K to L.}

assert K eq Domain(h); assert L eq Codomain(h);
if Type(f) eq RngUPolElt then
  coeffs := Coefficients(f);
  coeffs0 := [ CoerceToSubfieldElement(c, L, K, h) : c in coeffs ];
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
{Coerces the matrix M to the subfield K of L via the map hKL.}

return Matrix(K, [ [ CoerceToSubfieldElement(a, L, K, h) : a in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic CompositumExtra(K::Fld, L::Fld : CompatiblePlace:=false) -> Fld
{Returns the compositum of the fields K and L over their common base field, by taking a root of the polynomial defining L over K.}

// Trivial cases
if IsQQ(K) and IsQQ(L) then
  M := K;
  hKM := hom<K -> M | >;
  hLM := hom<L -> M | >;
  return M, hKM, hLM;
elif IsQQ(K) then
  M := L;
  hKM := hom<K -> M | >;
  hLM := hom<L -> M | M.1>;
  return M, hKM, hLM;
elif IsQQ(L) then
  M := K;
  hKM := hom<K -> M | M.1>;
  hLM := hom<L -> M | >;
  return M, hKM, hLM;
end if;

F, _, FinK, hFK := InclusionOfBaseExtra(K);
F, _, FinL, hFL := InclusionOfBaseExtra(L);

// Find minimal polynomial of both fields over the base
fK := MinimalPolynomial(K.1, FinK);
fL := MinimalPolynomial(L.1, FinL);
fLF := CoerceToSubfieldPolynomial(fL, FinL, F, hFL);
fKF := CoerceToSubfieldPolynomial(fK, FinK, F, hFK);

// Take corresponding root
fLK := ConjugatePolynomial(hFK, fLF);
g := Factorization(fLK, K)[1][1];
M, rL, hKM := NumberFieldExtra(g);
hLM := hom<L -> M | rL >;

if not CompatiblePlace then return M, hKM, hLM; end if;

// Take place compatible with both previous ones
CC := K`CC; genK := hKM(K.1); genL := hLM(L.1);
genKCC0 := EmbedExtra(K.1);
genLCC0 := EmbedExtra(L.1);
for iotaM in InfinitePlacesExtra(M) do
  genKCC := EmbedExtra(genK : iota := iotaM);
  genLCC := EmbedExtra(genL : iota := iotaM);
  if Abs(genKCC - genKCC0) lt CC`epscomp then
    if Abs(genLCC - genLCC0) lt CC`epscomp then
      M`iota := iotaM;
      return M, hKM, hLM;
    end if;
  end if;
end for;
error "No compatible infinite place found while taking compositum";

end intrinsic;


intrinsic CompositumExtra(Ks::List) -> Fld, List
{Returns the compositum of the fields in Ks over their common base field, together with inclusion maps.}

if #Ks eq 1 then
  return Ks[1], [* CanonicalInclusionMap(Ks[1], Ks[1]) *];
elif #Ks eq 2 then
  M, phi1, phi2 := CompositumExtra(Ks[1], Ks[2]);
  return M, [* phi1, phi2 *];
end if;
L, psis := CompositumExtra(Ks[1..(#Ks - 1)]);
M, phi1, phi2 := CompositumExtra(L, Ks[#Ks]);
return M, [* psi*phi1 : psi in psis *] cat [* phi2 *];

end intrinsic;
