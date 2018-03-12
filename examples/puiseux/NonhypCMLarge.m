AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

QQ := Rationals();
S<x,y,z> := PolynomialRing(QQ, 3);
fX := 3278898472*x^4 + 35774613556*x^3*y - 172165788624*x^3*z -
42633841878*x^2*y^2 + 224611458828*x^2*y*z + 362086824567*x^2*z^2 +
6739276447*x*y^3 + 195387780024*x*y^2*z + 1153791743988*x*y*z^2 -
3461357269578*x*z^3 - 18110161476*y^4 - 549025255626*y^3*z -
482663555556*y^2*z^2 + 15534718882176*y*z^3 - 61875497274721*z^4;
R<t> := PolynomialRing(QQ);
f := t^2 - t + 5;
h := hom< S -> R | [0, t, 1] >;
g := h(fX);
F<r> := Compositum(NumberField(f), NumberField(g));

//p := t;
//F<r> := NumberField(p);
rP := Roots(g, F)[1][1];
r2 := Roots(f, F)[1][1];

/*
QQ := Rationals();
S<x,y,z> := PolynomialRing(QQ, 3);
fX := 3278898472*x^4 + 35774613556*x^3*y - 172165788624*x^3*z -
42633841878*x^2*y^2 + 224611458828*x^2*y*z + 362086824567*x^2*z^2 +
6739276447*x*y^3 + 195387780024*x*y^2*z + 1153791743988*x*y*z^2 -
3461357269578*x*z^3 - 18110161476*y^4 - 549025255626*y^3*z -
482663555556*y^2*z^2 + 15534718882176*y*z^3 - 61875497274721*z^4;
QQ := Rationals();
R<t> := PolynomialRing(QQ);
f := t^6 - t^5 + 2*t^4 + 255*t^3 + 1291*t^2 - 784*t + 2192;
h := hom< S -> R | [0, t, 1] >;
g := h(fX);
F<r> := SplittingField(f*g);
rP := Roots(g, F)[1][1];
r6 := Roots(f, F)[1][1];
r2 := 1/19752*(25*r6^5 - 23*r6^4 + 114*r6^3 + 5199*r6^2 + 33481*r6 + 3028);
*/

P2<x,y,z> := ProjectiveSpace(F, 2);
fX := 3278898472*x^4 + 35774613556*x^3*y - 172165788624*x^3*z -
42633841878*x^2*y^2 + 224611458828*x^2*y*z + 362086824567*x^2*z^2 +
6739276447*x*y^3 + 195387780024*x*y^2*z + 1153791743988*x*y*z^2 -
3461357269578*x*z^3 - 18110161476*y^4 - 549025255626*y^3*z -
482663555556*y^2*z^2 + 15534718882176*y*z^3 - 61875497274721*z^4;
X := Curve(P2, fX);
P0 := X ! [0, rP, 1];

M := Matrix(F, [
[r2, 0, 0],
[0, r2, 0],
[0, 0, 1 - r2]
]);

/*
M := Matrix(F, [
[-r6, 1/1639416*(-1763*r6^5 + 2840*r6^4 + 9540*r6^3 - 529785*r6^2 - 589622*r6 + 3617300), 0],
[0, 1/819708*(-1763*r6^5 + 2840*r6^4 + 9540*r6^3 - 529785*r6^2 - 1409330*r6 + 3617300), 0],
[0, 0, 1/1639416*(-4329*r6^5 - 4379*r6^4 - 7494*r6^3 - 1198119*r6^2 - 7501575*r6 - 5997476)]
]);

M := Matrix(F, [
[r2, 0, 0],
[0, r2, 0],
[0, 0, 1 - r2]
]);
*/

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Tangent representation:";
print M;

print "Calculating divisor...";
//time test, D := DivisorFromMatrixRRSplit(X, P0, X, P0, M : LowerBound := 168, Margin := 2^8);
time test, D := DivisorFromMatrixRRSplit(X, P0, X, P0, M : LowerBound := 30, Margin := 2^8);
print D;

exit;

print "Irreducible components:";
print IrreducibleComponents(D);

exit;
