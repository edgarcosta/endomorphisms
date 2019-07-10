//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 2);

prec := 300;

F := RationalsExtra(prec);
P3<x,y,z,w> := ProjectiveSpace(F, 3);

f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w
+ 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1, f2]);

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^6 - t^3 + 1, prec);
S<t> := PolynomialRing(F);
g := t^3 + (1212597*r^4 + 989496*r)*t - 259785954*r^3 + 1111959306;
Y := HyperellipticCurve(g);

print "";
print "Curve X:";
print X;

print "";
print "Curve Y:";
print Y;

P := PeriodMatrix(X);
Q := PeriodMatrix(Y);
GeoHomRep := GeometricHomomorphismRepresentation(P, Q, F);

print "";
print "Homomorphisms:";
print GeoHomRep;

print "";
print "Degrees:";
EX := StandardSymplecticMatrix(4);
EY := StandardSymplecticMatrix(1);
Rs := [ tup[2] : tup in GeoHomRep ];
for R in Rs do
    test, d := IsRationalMultiple(R*EX*Transpose(R), EY);
    assert test;
    print d;
end for;

