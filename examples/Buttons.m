SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 100;

// This one takes quite some time!
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^6 + r;

// More examples over an extension
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^5 + x + 1;
f := x^5 + r*x^3 + x;
X := HyperellipticCurve(f);

// More examples over an extension
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
R<x> := PolynomialRing(F);
f := x^6 + r;
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ];
X := HyperellipticCurve(f);

// More examples over QQ
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
// Big Sato-Tate group, this calculation takes about 20 minutes:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8;
// Decomposition:
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x);
X := HyperellipticCurve(f);

// Squares: only enable if you have curve_reconstruction
// TO BE DEBUGGED f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2;
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10;
X := HyperellipticCurve(f);

// Plane curves
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);

f := y^3*z - x^4 - z^4;
f := 1 + 7*x*y + 21*x^2*y^2 + 35*x^3*y^3 + 28*x^4*y^4 + 2*x^7 + 2*y^7;
f := x^3*z + 2*x^2*y^2 + x^2*y*z + 2*x^2*z^2 - x*y^2*z + x*y*z^2 - x*z^3 + y^3*z - y^2*z^2 + y*z^3 - z^4;
f := x^3*z + x^2*y^2 - 3*x*y^2*z - 4*x*z^3 - 2*y^4 + y^3*z - 4*y^2*z^2 - 3*z^4;
f := x^3*z + 2*x^2*y^2 + x^2*y*z + 3*x^2*z^2 - 4*x*y^3 - x*y^2*z + 5*x*y*z^2 + x*z^3 + 2*y^4 + 6*y^3*z + 6*y^2*z^2 + 2*y*z^3;
f := 2*x^4 + 3*x^3*y + 4*x^3*z + 6*x^2*y^2 + 4*x^2*y*z + 7*x^2*z^2 + 4*x*y^3 + 4*x*y^2*z + 7*x*y*z^2 + 4*x*z^3 + 3*y^4 + 2*y^3*z + 3*y^2*z^2 + 5*y*z^3 + 2*z^4;
f := x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + 4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4;
X := PlaneCurve(f);

// General curve
F := RationalsExtra(prec);
P3<x,y,z,w> :=ProjectiveSpace(F, 3);

f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w + 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1, f2]);

// More examples over QQ
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := -25*x^6 + 12*x^5 + 27*x^4 - 16*x^3 - 3*x^2 + 4*x + 1;
X := HyperellipticCurve(f);


print "";
print "Curve:";
print X;

desc := HeuristicEndomorphismDescription(X);
print "";
print "Heuristic endomorphism description:";
print desc;

A := HeuristicEndomorphismAlgebra(X);
print "";
print "Heuristic endomorphism algebra:";
print A;

OO := HeuristicEndomorphismRing(X);
print "";
print "Heuristic endomorphism ring:";
print OO;

rep := HeuristicEndomorphismRepresentation(X);
print "";
print "Heuristic endomorphism representation:";
print rep;

L := HeuristicEndomorphismFieldOfDefinition(X);
print "";
print "Heuristic field of definition of the endomorphisms:";
print L;

Lat := HeuristicEndomorphismLattice(X);
print "";
print "Heuristic endomorphism lattice:";
print Lat;

test_gl2 := HeuristicIsGL2(X);
print "";
print "Heuristic GL_2-determination:";
print test_gl2;

/*
facinfo := HeuristicJacobianFactors(X);
print "";
print "Heuristic Jacobian factors:";
print facinfo;

exps, test, degs := IsogenyInformation(X : facinfo := facinfo);
print "";
print "Isogeny exponents:";
print exps;
print "";
print "Compatible with polarizations:";
print test;
print "";
print "Degrees (if applicable):";
print degs;
*/

