/*
  An example in Magma.
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 300;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

// Big Sato-Tate group, this calculation takes about 20 minutes:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8;

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^5 + x + 1;
f := x^5 + r*x^3 + x;

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
//F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ];
f := x^6 + r;

X := HyperellipticCurve(f);

F := RationalsExtra(prec);
P3<x,y,z,w> :=ProjectiveSpace(F, 3);

f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w
+ 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1,f2]);

F := RationalsExtra(prec); R<x> := PolynomialRing(F); X := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]);

print "";
print "Curve:";
print X;

print "";
print "Heuristic field of definition of the endomorphisms:";
L := HeuristicEndomorphismFieldOfDefinition(X);
print L;

print "";
print "Heuristic geometric endomorphism algebra:";
A := HeuristicEndomorphismAlgebra(X : Geometric := true);
print A;
print "Description:";
desc := HeuristicEndomorphismAlgebraDescription(X : Geometric := true);
print desc;

print "";
print "Heuristic endomorphism algebra over the base:";
A := HeuristicEndomorphismAlgebra(X);
print A;
print "Description:";
desc := HeuristicEndomorphismAlgebraDescription(X);
print desc;

print "";
print "Heuristic endomorphism lattice:";
Lat := HeuristicEndomorphismLattice(X);
print Lat;

print "";
print "Heuristic GL_2-determination (after Ribet):";
test_gl2_ribet := HeuristicIsGL2Ribet(X);
print test_gl2_ribet;

print "";
print "Heuristic GL_2-determination (generalized notion):";
test_gl2_gen := HeuristicIsGL2Generalized(X);
print test_gl2_gen;

print "";
print "Heuristic Jacobian factors:";
facs := HeuristicJacobianFactors(X);
print facs;
