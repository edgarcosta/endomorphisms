AttachSpec("../../endomorphisms/magma/spec");

print "Over QQ:";

F := Rationals();
R<x> := PolynomialRing(F);
f := x^3 + x + 1;
g := x^2 - 5;

K := NumberFieldExtra(f);
SetInfinitePlace(K, InfinitePlaces(K)[2]);

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
//print L`iota;
//print Roots(f, L);


print "Over a number field:";

F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);
f := x^3 + x - r;
g := x^2 - 5;
g := x^2 - 5 + r;

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
print Roots(f, L);

print L;
print MinimalPolynomial(L.1, F);

Gp, Gf, Gphi := AutomorphismGroup(L);
print RelativeFixedField(L, [ Gphi(Gp.1) ]);
print RelativeFixedField(L, [ Gphi(g) : g in Gp ]);

K := RelativeFixedField(L, [ Gphi(Gp.2) ]);
print K;
print K ! L ! K.1;

pol := MinimalPolynomial(K.1/14, F);
K := NumberField(pol);
K := ClearFieldDenominator(K);
print K;
print MinimalPolynomial(K.1, F);
print MinimalPolynomial(K.1);

/*
SetInfinitePlace(L, InfinitePlaces(L)[7]);
print L`iota;
print F`iota;
*/

print "Linear extension:";

F := Rationals();
R<x> := PolynomialRing(F);

print RelativeSplittingFieldExtra(x);
print RelativeSplittingFieldExtra(x^2 - 2);

F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);

L := RelativeSplittingFieldExtra(x);
print L;
Gp, Gf, Gphi := AutomorphismGroup(L);
print Gp;

exit;
