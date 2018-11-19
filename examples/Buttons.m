/*
  An example in Magma.
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 500;
// This one takes quite some time!
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^6 + r;


F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
// Big Sato-Tate group, this calculation takes about 20 minutes:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8;
// Squares:
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10; h := R ! 0;
X := HyperellipticCurve(f);


/*
F := RationalsExtra(prec);
R<x,y> := PolynomialRing(F, 2);
z := 1;

f := y^3*z - x^4 - z^4;
f := x^3*z + 2*x^2*y^2 + x^2*y*z + 2*x^2*z^2 - x*y^2*z + x*y*z^2 - x*z^3 + y^3*z - y^2*z^2 + y*z^3 - z^4;
f := x^3*z + x^2*y^2 - 3*x*y^2*z - 4*x*z^3 - 2*y^4 + y^3*z - 4*y^2*z^2 - 3*z^4;
f := x^3*z + 2*x^2*y^2 + x^2*y*z + 3*x^2*z^2 - 4*x*y^3 - x*y^2*z + 5*x*y*z^2 + x*z^3 + 2*y^4 + 6*y^3*z + 6*y^2*z^2 + 2*y*z^3;
f := 2*x^4 + 3*x^3*y + 4*x^3*z + 6*x^2*y^2 + 4*x^2*y*z + 7*x^2*z^2 + 4*x*y^3 + 4*x*y^2*z + 7*x*y*z^2 + 4*x*z^3 + 3*y^4 + 2*y^3*z + 3*y^2*z^2 + 5*y*z^3 + 2*z^4;
X := PlaneCurve(f);


F := RationalsExtra(prec);
P3<x,y,z,w> :=ProjectiveSpace(F, 3);

f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w
+ 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1, f2]);
*/


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
facs := HeuristicJacobianFactors(X : AllMaps := false);
print facs;
