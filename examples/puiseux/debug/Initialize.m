import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;

F := RationalsExtra();
R<x> := PolynomialRing(F);
X := HyperellipticCurve(x^6 + 3*x + 4);
P := X ! [ 0, 2, 1 ];
P := X ! [ 1, 1, 0 ];
//X := HyperellipticCurve(6*x^5 - 6*x^4 - 2*x^3 + 5*x^2 - 1);
//P := X ! [ 1, 0, 0 ];
InitializeCurve(X, P : AssertNonWP := true);
print "";
print X`U;
print X`DEs;
print X`unif_index;
print X`patch_index;

S<x,y,z> := PolynomialRing(F, 3);
X := _PlaneCurve(3*x^3*y + 5*y^3*z + 7*z^3*x);
P := X ! [ 0, 1, 0 ];
InitializeCurve(X, P);
print "";
print X`U;
print X`DEs;
print X`unif_index;
print X`patch_index;

S<x,y,z> := PolynomialRing(F, 3);
X := _PlaneCurve(y^4 + y^2*(5*x^2 - x*z - 3*z^2) + x*z*(3*x^2 + 2*x*z - 7*z^2));
X := _PlaneCurve(x^4 + x^2*(5*y^2 - y*z - 3*z^2) + y*z*(3*y^2 + 2*y*z - 7*z^2));
P := X ! [ 0, 0, 1 ];
InitializeCurve(X, P);
print "";
print X`U;
print X`DEs;
print X`unif_index;
print X`patch_index;

/* NOTE: This means that putting x first is somewhat bad in this case. However,
 * this should be done by hand; Magma's differential functionality puts d before
 * the first coordinate, and if we do not take a uniformizer there, we lose
 * precision. */
print "";
K<u,v> := X`KU;
print v^4;

