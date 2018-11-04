/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

//SetVerbose("EndoFind", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

// Big Sato-Tate group, this calculation takes about 20 minutes:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1; h := R ! 0;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8; h := 0;

R<t> := PolynomialRing(Rationals());
F<r> := NumberFieldExtra(t^2 - 5);
R<x> := PolynomialRing(F);
f := x^5 + x + 1; h := R ! 0;
f := x^5 + r*x^3 + x; h := R ! 0;

R<t> := PolynomialRing(Rationals());
F<r> := NumberFieldExtra(t^2 - t + 1);
R<x> := PolynomialRing(F);
f := x^6 + r; h := R ! 0;
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ]; h := R ! 0;

X := HyperellipticCurve(f, h);
print "Curve:";
print X;

P := PeriodMatrix(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

/* Entries can be made relative by using RelativeField if so desired */
print "";
print "Geometric endomorphism representations:";
print GeoEndoRep;

Rs := [ tup[2] : tup in GeoEndoRep ];
Rs := &cat[ [ BlockMatrix([[R,0],[0,0]]), BlockMatrix([[0,R],[0,0]]), BlockMatrix([[0,0],[R,0]]), BlockMatrix([[0,0],[0,R]]) ] : R in Rs ];

A := Algebra(MatrixRing(Rationals(), 2*4));
GensA := [ A ! Eltseq(R) : R in Rs ];
/* As a subalgebra */
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
/* As an associative algebra */
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

/* Make heuristic algorithm using John's remarks on characteristic polynomials:
 * eliminate primes outside discriminant of order an use heuristic methods there
 * */
/* Implement isomorphism test in linked article */

/* Find symmetric idempotents quickly */

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f := x^4 + x^2; h := x^3 + 1;
f := x^6 + 3*x^5 + 6*x^4 + 7*x^3 + 6*x^2 + 3*x + 1; h := x^2 + x;

X := HyperellipticCurve(f, h);
print "Curve:";
print X;

P := PeriodMatrix(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

/* Entries can be made relative by using RelativeField if so desired */
print "";
print "Geometric endomorphism representations:";
print GeoEndoRep;

GeoEndoRepFixed := RosatiFixedModule(GeoEndoRep);
Rs := [ tup[2] : tup in GeoEndoRepFixed ];

S := PolynomialRing(Rationals(), #Rs);
Rs := [ ChangeRing(R, S) : R in Rs ];
M := &+[ S.i*Rs[i] : i in [1..#Rs] ];
D := M^2 - M;
S := Scheme(AffineSpace(S), Eltseq(D));

print "";
print "Symmetric idempotents";
print IrreducibleComponents(S);
