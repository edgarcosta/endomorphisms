AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

// Example 8
// Interesting because of degrees involved
R<t> := PolynomialRing(Rationals());
f := x^6 + 42*x^4 + 441*x^2 + 784;
F<s> := NumberField(f);
R<x> := PolynomialRing(F);
//R<x> := PolynomialRing(Rationals());
f := 16*x^7 + 357*x^5 - 819*x^3 + 448*x;
f /:= 2;
X := HyperellipticCurve(f);
//print RationalPoints(X : Bound := 100);
print [ Evaluate(f, n) : n in [-3..3] ];
P0 := X ! [1, 1, 1];

Ms := [ ];

//M := Matrix(F, [
//[1/28*(-s^3 - 21*s),                  0,                  0],
//[                 0,  1/28*(s^3 + 21*s),                  0],
//[                 0,                  0, 1/28*(-s^3 - 21*s)]
//]);
//Append(~Ms, M);

M := Matrix(F, [
[1/168*(-5*s^4 - 133*s^2 - 392), 0, 1/168*(s^4 - 7*s^2 - 392)],
[0, 1/56*(3*s^4 + 91*s^2 + 392), 0],
[1/84*(s^4 - 7*s^2 - 392), 0, 1/42*(-s^4 - 35*s^2 - 196)]
]);
Append(~Ms, M);

//M := Matrix(F, [
//[1/168*(-s^4 + 7*s^2 + 392), 0, 1/42*(-s^4 - 35*s^2 - 196)],
//[0, 1/28*(s^4 + 21*s^2), 0],
//[1/21*(-s^4 - 35*s^2 - 196), 0, 1/168*(-5*s^4 - 133*s^2 - 392)]
//]);
//Append(~Ms, M);

//M := Matrix(F, [
//[1/168*(s^5 + 35*s^3 + 154*s), 0, 1/168*(s^5 + 35*s^3 + 322*s)],
//[0, 1/56*(s^5 + 35*s^3 + 210*s), 0],
//[1/84*(s^5 + 35*s^3 + 322*s), 0, 1/84*(s^5 + 35*s^3 + 238*s)]
//]);
//Append(~Ms, M);

print "Field:";
print F;
print "Curve:";
print X;

for M in Ms do
    print "Tangent representation:";
    print M;
    print "Minimal polynomial:";
    print MinimalPolynomial(M);

    print "Calculating divisor:";
    time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
    eqs := DefiningEquations(D);
    R<y2,y1,x2,x1> := Parent(eqs[1]);
    print "Divisor:";
    print D;
end for;

exit;


AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

// Weng 2
// Verified in a few minutes
R<t> := PolynomialRing(Rationals());
f := x^6 + 6*x^4 + 9*x^2 + 1;
F<s> := NumberField(f);
R<x> := PolynomialRing(F);
f := x^7 + 6*x^5 + 9*x^3 + x;
f /:= 17;
X := HyperellipticCurve(f);
P0 := X ! [1, 1];

Ms := [ ];

M := Matrix(F, [
[    1/3*(s^4 + 2*s^2 - 2),                         0,  1/3*(-2*s^4 - 7*s^2 - 2)],
[                        0,                   s^2 + 2,                         0],
[1/3*(-4*s^4 - 14*s^2 - 4),                         0,    1/3*(-s^4 - 5*s^2 - 4)]
]);
Append(~Ms, M);

M := Matrix(F, [
[ s^3 + 3*s,          0,          0],
[         0, -s^3 - 3*s,          0],
[         0,          0,  s^3 + 3*s]
]);
Append(~Ms, M);

M := Matrix(F, [
[1/3*(-s^5 - 5*s^3 - 4*s), 0, 1/3*(2*s^5 + 10*s^3 + 11*s)],
[0, s, 0],
[1/3*(4*s^5 + 20*s^3 + 22*s), 0, 1/3*(s^5 + 5*s^3 + 7*s)]
]);
Append(~Ms, M);

print "Field:";
print F;
print "Curve:";
print X;

for M in Ms do
    print "Tangent representation:";
    print M;
    print "Minimal polynomial:";
    print MinimalPolynomial(M);

    print "Calculating divisor:";
    time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
    eqs := DefiningEquations(D);
    R<y2,y1,x2,x1> := Parent(eqs[1]);
    print "Divisor:";
    print D;
end for;

exit;


AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

// Weng 3
// Verified; it takes about 6 hours and 6 Gb of memory
R<t> := PolynomialRing(Rationals());
f := x^6 + 13*x^4 + 50*x^2 + 49;
F<r> := NumberField(f);
R<x> := PolynomialRing(F);
f := 21*x^7 + 37506*x^5 + 933261*x^3 + 5841759*x;
X := HyperellipticCurve(f);
P0 := X ! [3, -7203];

M := Matrix(F, [
[1/7*(r^5 + 6*r^3 + r), 0, 0],
[0, 1/7*(-r^5 - 6*r^3 - r), 0],
[0, 0, 1/7*(r^5 + 6*r^3 + r)]
]);
M := Transpose(M);

M := Matrix(F, [
[1/49*(r^5 - r^3 - 13*r), 0, 1/49*(660*r^5 + 5808*r^3 + 10824*r)],
[0, r, 0],
[1/539*(5*r^5 + 44*r^3 + 82*r), 0, 1/49*(6*r^5 + 43*r^3 + 69*r)]
]);
M := Transpose(M);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;
print "Minimal polynomial:";
print MinimalPolynomial(M);

print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;


AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

// Not verified yet
R<t> := PolynomialRing(Rationals());
f := x^6 - 14*x^3 + 63*x^2 + 168*x + 161;
F<s> := NumberField(f);
R<x> := PolynomialRing(F);
f := -4*x^8 + 105*x^6 - 945*x^4 + 2100*x^2 - 5895*x + 420;
h := x^4;
p := (4*f + h^2)/1680;
X := HyperellipticCurve(p);
P0 := X ! [0, 1];

M := Matrix(F, [
[1/506*(-9*s^5 + 12*s^4 - 16*s^3 + 63*s^2 - 651*s - 897), 0, 0],
[0, 1/506*(9*s^5 - 12*s^4 + 16*s^3 - 63*s^2 + 651*s + 391), 0],
[1/253*(-18*s^5 + 24*s^4 - 32*s^3 + 126*s^2 - 1302*s - 1288), 0, 1/506*(9*s^5 - 12*s^4 + 16*s^3 - 63*s^2 + 651*s + 391)]
]);

M := Matrix(F, [
[1/253*(-8*s^5 + 3*s^4 + 19*s^3 + 102*s^2 - 663*s - 1127), 0, 0],
[1/3795*(-4*s^5 + 82*s^4 - 2*s^3 + 74*s^2 - 964*s + 2898), 1/2530*(-22*s^5 - 9*s^4 - 126*s^3 + 637*s^2 - 1507*s - 2576), 1/3795*(2*s^5 - 41*s^4 + s^3 - 37*s^2 + 482*s - 1449)],
[1/1265*(-52*s^5 - 84*s^4 + 319*s^3 + 272*s^2 - 3677*s - 13041), 1/1265*(6*s^5 - 123*s^4 + 3*s^3 - 111*s^2 + 1446*s - 4347), 1/2530*(-28*s^5 + 114*s^4 - 129*s^3 + 748*s^2 - 2953*s + 1771)]
]);

M := Matrix(F, [
[1/253*(-8*s^5 + 3*s^4 + 19*s^3 + 102*s^2 - 410*s - 1127), 0, 0],
[1/3795*(28*s^5 - 22*s^4 + 152*s^3 - 794*s^2 + 2194*s + 1932), 1/2530*(-36*s^5 + 94*s^4 - 179*s^3 + 988*s^2 - 3363*s + 161), 1/3795*(-14*s^5 + 11*s^4 - 76*s^3 + 397*s^2 - 1097*s - 966)],
[1/1265*(-86*s^5 - 31*s^4 + 141*s^3 + 1223*s^2 - 4028*s - 14329), 1/1265*(-42*s^5 + 33*s^4 - 228*s^3 + 1191*s^2 - 3291*s - 2898), 1/2530*(6*s^5 + 61*s^4 + 49*s^3 - 203*s^2 - 72*s + 3059)]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
