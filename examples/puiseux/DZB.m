AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^6 - t^3 + 1);
P2<x,y,z> := ProjectiveSpace(F, 2);

fX := x^3*z^2 + 108*x^2*y^3 + 36*x^2*y*z^2 + 96*x^2*z^3 + 2592*x*y^4 + 8208*x*y^3*z + 432*x*y^2*z^2 + 2304*x*y*z^3 + 2928*x*z^4 + 15552*y^5 + 98496*y^4*z + 157248*y^3*z^2 + 13824*y^2*z^3 + 35136*y*z^4 + 27584*z^5;
X := Curve(P2, fX);
P0 := X ! [ -4*r^5 - 4*r^4 + 8*r^2 + 8*r - 32, 0, 1 ];

S<t> := PolynomialRing(F);
fY := t^3 + (1212597*r^4 + 989496*r)*t - 259785954*r^3 + 1111959306;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [ [1/91*(-38*r^3 + 22), 1/91*(-18*r^4 + 20*r), 1/91*(-216*r^4 - 1456*r^3 + 240*r + 728), 1/91*(-640*r^4 + 792*r)] ]);
M := Matrix(F, [ [1/91*(-22*r^3 - 16), 1/91*(2*r^4 + 18*r), 1/91*(24*r^4 - 728*r^3 + 216*r - 728), 1/91*(152*r^4 + 640*r)] ]);
M := Matrix(F, [ [1/91*(-54*r^4 + 60*r), 1/91*(16*r^3 - 38), 1/91*(-1920*r^4 + 192*r^3 + 2376*r - 456), 8*r^3 - 16] ]);
M := Matrix(F, [ [1/91*(-6*r^4 + 38*r^3 - 54*r - 22), 1/91*(18*r^4 - 38*r^3 - 20*r + 22), 1/91*(-240*r^4 + 1000*r^3 - 2160*r - 464), 1/91*(640*r^4 - 1456*r^3 - 792*r + 728)] ]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Curve Y:";
print Y;
print "Point Q0:";
print Q0;
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
print "Original equation:", fX;
print "First affine patch:", R ! DefiningEquation(AffinePatch(X, 1));

/* Check that the answer is a projection after transformation, automatic: */
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;
