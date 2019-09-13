SetVerbose("EndoFind", 3);
SetVerbose("CurveRec", 0);

prec := 100;
F := RationalsExtra(prec);
CC := F`CC;
R<x> := PolynomialRing(F);

f3 := x^8 + x^6 + 4*x^4 - 7*x^2 - 17;
f1 := x^4 + x^3 + 4*x^2 - 7*x^1 - 17;
f2 := x*f1;

X3 := HyperellipticCurve(f3);
X2 := HyperellipticCurve(f2);
X1 := HyperellipticCurve(f1);

P3 := PeriodMatrix(X3);
P2 := PeriodMatrix(X2);
P1 := PeriodMatrix(X1);

function RandomInvertibleMatrix(g, B)
D := [ -B..B ];
repeat
    T := Matrix(Rationals(), g,g, [ Random(D) : i in [1..g^2] ]);
until Determinant(T) eq 1;
return T;
end function;

U3 := RandomInvertibleMatrix(6, 2);
U3CC := ChangeRing(U3, BaseRing(P3));
Q3 := P3*U3CC;
print "";
print U3;

U2 := RandomInvertibleMatrix(4, 2);
U2CC := ChangeRing(U2, BaseRing(P2));
Q2 := P2*U2CC;
print "";
print U2;

test, E3 := SomePrincipalPolarization(P3);
print E3;

/* Test projection */
GeoHomRep := GeometricHomomorphismRepresentationCC(P3, P2);
A, R := Explode(GeoHomRep[1]);
print "";
print R;

test, E3 := SomePrincipalPolarization(Q3);
_, U := FrobeniusFormAlternatingRational(E3);
print SomePrincipalPolarization(Q3*ChangeRing(Transpose(U), BaseRing(Q3)));

print "";
print IsBigPeriodMatrix(P3);

print "";
print SymplecticSubmodules(12, 1);

E := Matrix(Rationals(), [[0,1],[-1,0]]);
E := Matrix(Rationals(), [[0,0,1,0],[0,0,0,12],[-1,0,0,0],[0,-12,0,0]]);
U2 := RandomInvertibleMatrix(4, 2);
E := U2*E*Transpose(U2);
print "";
Ts := IsogenousPPLattices(E : ProjToPP := true);

for T in Ts do
    print "";
    print T^(-1);
    print Determinant(T^(-1));
    print T*E*Transpose(T);
end for;



