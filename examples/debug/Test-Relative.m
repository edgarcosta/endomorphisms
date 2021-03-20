SetVerbose("EndoFind", 0);

prec := 200;
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
K<r> := NumberFieldExtra(x^2 - 2);
S<y> := PolynomialRing(K);
L<s> := NumberFieldExtra(y^2 - 3);
Lalt := NumberFieldExtra(y^2 - 27);

print L;
print L`base;
print L`base_gen;
print L`CC;
print L`iota;

print "";
print FieldDescriptionExtra(L);
print ElementDescriptionExtra(L.1);
print InfinitePlacesExtra(L);
print EmbedExtra(r);

print "";
print EmbedPolynomialsExtra([ y^2 - 3, y^2 - 5 ]);

print "";
Gp, Gf, Gphi := AutomorphismGroupPari(L);
gens := [ Gphi(sigma) : sigma in Gp ];
K, h := FixedFieldExtra(L, gens);
print K;
print K`base;
print K`base_gen;
print K`CC;
print K`iota;
print h;

print "";
test, h := IsIsomorphic(L, Lalt);
TransferAttributesExtra(L, Lalt, h);
print Lalt;
print Lalt`base;
print Lalt`base_gen;
print Lalt`CC;
print Lalt`iota;
print h;

print "";
aCC := Sqrt(CC ! 2);
print MinimalPolynomialExtra(aCC, K);

aCC := Sqrt(CC ! 3);
print MinimalPolynomialExtra(aCC, K);

print "";
time K, elts := NumberFieldExtra([ Sqrt(CC ! i ) : i in [1..5] ], F);
time K, elts := NumberFieldExtra([ Sqrt(CC ! i ) : i in [1..5] ], L);
print "";
print K;
print K`base;
print K`base_gen;
print K`CC;
print K`iota;
print elts;

print "";
time K, elts := SplittingFieldExtra([ Sqrt(CC ! i ) : i in [1..5] ], F);
time K, elts := SplittingFieldExtra([ Sqrt(CC ! i ) : i in [1..5] ], L);
print "";
print K;
print K`base;
print K`base_gen;
print K`CC;
print K`iota;
print elts;

h := hom< F -> K | >;
print CoerceToSubfieldElement(K`base_gen, K, F, h);
print CoerceToSubfieldElement(K ! 1, K, F, h);

F := RationalsExtra(prec);
h := hom< F -> K | >;
print CoerceToSubfieldElement(K ! 1, K, F, h);

R<x> := PolynomialRing(Rationals());
K<r> := BaseNumberFieldExtra(x^2 - 2, prec);
S<y> := PolynomialRing(K);
L1, s := NumberFieldExtra(y^2 - (r + 1));
L2, s := NumberFieldExtra(y^2 - (r + 2));
L3, s := NumberFieldExtra(y^3 + y - (r + 3));

M := CompositumExtra(L1, L2);
print M;
M := CompositumExtra([* L1, L2, L3 *]);
print M;
M := CompositumExtra([* L1, L1 *]);
print M;


// We work at precision 300
prec := 300;

// Creation of base number field (say the CM field of its reflex field).
// All other number fields henceforth will be over this base field.
R<t> := PolynomialRing(Rationals());
f := x^3 - 2;
K := BaseNumberFieldExtra(f, prec);

print "";
print "Number field:";
print K;
print "";
print "Embedding of K into CC defined by:";
print K`iota;

// Associated complex field
CC := K`CC;
// We create some algebraic numbers in CC
a1 := Sqrt(CC ! 2);
a2 := Sqrt(CC ! 3);
a3 := Sqrt(CC ! 5);

// Adjoining the elements and checking that they embed correctly
L1, as1 := NumberFieldExtra([ a1, a2 ], K);
print "";
print "Number field after first adjunction:";
print L1;
print "Embedding of L1 into CC defined by:";
print L1`iota;
print "Algebraizations of given complex numbers:";
print as1;
print "Embeddings of algebraizations:";
print [ EmbedExtra(a) : a in as1 ];


// This can be repeated, in the sense that a3 can still be adjoined afterwards
L2, as2, h := NumberFieldExtra([ a3 ], L1);
print "";
print "Number field after second adjunction:";
print L2;
print "Embedding of L2 into CC defined by:";
print L2`iota;
print "Algebraizations of given complex numbers:";
print as2;
print "Embeddings of algebraizations:";
print [ EmbedExtra(a) : a in as2 ];
print "Embeddings of previous algebraizations:";
print [ EmbedExtra(h(a)) : a in as1 ];



