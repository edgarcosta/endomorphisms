SetVerbose("EndoFind", 1);

R<x> := PolynomialRing(Rationals());
fs := [ x^2 + d : d in [1..3] ];
L := RelativeSplittingFieldQQ(fs);
print L;
print FieldDescription(L);
print ElementDescription(L.1);
SetInfinitePlaceUpwards(L);
print L`iota;

print "";
K<r> := NumberField(x^2 + 1);
S<y> := PolynomialRing(K);
L<s> := NumberField(y^2 + 2);
M := ExtendRelativeNumberField(y^2 + 3);
print M;

print "";
T<z> := PolynomialRing(L);
M := ExtendRelativeNumberField(z^2 + 3);
print M;

print "";
M := ExtendRelativeSplittingField(L, y^3 + 3);
print M;

print "";
M := RelativeSplittingFieldExtra([ x^2 + 2, x^2 + 3 ]);
print M;

print "";
M := RelativeSplittingFieldExtra([ y^2 + 2, y^2 + 3, y^2 + 5 ]);
print M;
print BaseRing(M);

/*
print "";
Ks := [* NumberField(x^2 + d) : d in [1..3] *];
print RelativeCompositum(Ks);
*/

print "";
Ls := [* NumberField(y^2 + d) : d in [2,3,5] *];
print RelativeCompositum(Ls[1], Ls[2]);
