SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 100;

Xs := [* *];

// More examples over an extension
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^5 + r*x^3 + x;
X := HyperellipticCurve(f); Append(~Xs, X);

// More examples over an extension
R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
R<x> := PolynomialRing(F);
f := x^6 + r;
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ];
X := HyperellipticCurve(f); Append(~Xs, X);

// More examples over QQ
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8;
X := HyperellipticCurve(f); Append(~Xs, X);

// Plane curves
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
fs := [ ];
f := x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + 4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4; Append(~fs, f);
for f in fs do
    X := _PlaneCurve(f); Append(~Xs, X);
end for;

// General curve
F := RationalsExtra(prec);
P3<x,y,z,w> := ProjectiveSpace(F, 3);
f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w + 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1, f2]); Append(~Xs, X);


for X in Xs[1..#Xs] do
    print "";
    print "Curve:";
    print X;

    time P := PeriodMatrix(X);
    //time desc := HeuristicEndomorphismAlgebra(X : CC := true);
    time desc := HeuristicEndomorphismAlgebra(X : Geometric := true);
    time rep := HeuristicEndomorphismRepresentation(X);
    time lat := HeuristicEndomorphismLattice(X);
    print lat;
end for;
