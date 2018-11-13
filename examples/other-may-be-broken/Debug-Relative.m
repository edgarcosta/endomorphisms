SetVerbose("EndoFind", 2);

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
print IsQQ(L);
print FieldDescription(L);
print ElementDescription(L.1);
print InfinitePlacesExtra(L);
print EvaluateExtra(r, L`iota);

print "";
print EmbedAtInfinitePlacePolynomials([ y^2 - 3, y^2 - 5 ]);


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

F := K`base;
h := hom< F -> K | K`base_gen >;
print CoerceToSubfieldElement(K`base_gen, K, F, h);
print CoerceToSubfieldElement(K ! 1, K, F, h);

F := RationalsExtra(prec);
h := hom< F -> K | >;
print CoerceToSubfieldElement(K ! 1, K, F, h);

prec := 200;
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
K<r> := NumberFieldExtra(x^2 - 2);
S<y> := PolynomialRing(K);
L1, s := NumberFieldExtra(y^2 - (r + 1));
L2, s := NumberFieldExtra(y^2 - (r + 2));
L3, s := NumberFieldExtra(y^3 + y - (r + 3));

print CompositumExtra(L1, L2);
print CompositumExtra([* L1, L2 *]);
M := CompositumExtra([* L1, L2, L3 *]);
print M;
