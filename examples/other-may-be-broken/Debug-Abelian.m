/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

prec := 600;
prec := 200;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := x^8 + x^6 + 2; h := R ! 0;
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
f := x^7 + x^6 + 5*x^5 - 3*x^4 + 2*x^3 - 13*x^2 + 7*x - 1; h := x^3 + x;
//f := x^7 + x^6 + x^5 + x^3 + x^2 + x; h := x^4 + x^2 + 1;
//f := -2*x^7 - 4*x^6 + 3*x^4 + x^3 - 2*x^2 - x; h := x^2 + x + 1;
//f := x^7 - 2*x^5 - 4*x^4 - 2*x^3 + x; h := x^4 + x^2 + 1;
//f := x^7 - 4*x^6 - x^5 + 10*x^4 + 3*x^3 - 6*x^2 - x + 1; h := x^3 + x;
//f := x^7 + x^6 - 4*x^5 + 4*x^3 - 5*x^2 + 2*x - 1; h := x^4 + x^3 + x + 1; SetEpsComp(F`CC, 10^70*F`CC`epscomp);

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
Y := HyperellipticCurve(x^3 - 7371/16*x - 120285/32);

P := PeriodMatrix(X);
Q := PeriodMatrix(Y);
HomRepPQ := GeometricHomomorphismRepresentationCC(P, Q);
HomRepQP := GeometricHomomorphismRepresentationCC(Q, P);
hPQ := HomRepPQ[1];
hQP := HomRepQP[1];

print hPQ[2];
print KerModKer0(hPQ, P, Q);
K, hK := Ker0(hPQ, P, Q);
print hK[2];
C, hC := Coker(hPQ, P, Q);
print hC[2];
I, hI := Img(hPQ, P, Q);
print hI[2];

print hQP[2];
print KerModKer0(hQP, Q, P);
K, hK := Ker0(hQP, Q, P);
print hK[2];
C, hC := Coker(hQP, Q, P);
print hC[2];
I, hI := Img(hQP, Q, P);
print hI[2];
