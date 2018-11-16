SetVerbose("EndoFind", 0);

print "Test SymplecticSubmodules:";
print #SymplecticSubmodules(4, 2);
print #SymplecticSubmodules(2, 2);
print #SymplecticSubmodules(2, 4);
//print #SymplecticSubmodules(2, 6);

print "Checking IsogenousPPLattices...";
d := 6;
EQ := Matrix(QQ, [[0,0,d,0],[0,0,0,1],[-d,0,0,0],[0,-1,0,0]]);
M := RandomSymplecticMatrix(2, 15);
EQ := M*EQ*Transpose(M);
Us := IsogenousPPLattices(EQ);
for U in Us do
    if not U*EQ*Transpose(U) eq d*ChangeRing(StandardSymplecticMatrix(2), Rationals()) then
        print "Problem";
        print U*EQ*Transpose(U);
    end if;
end for;
print "done.";

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f := x^8 + x^6 + 1;
X := HyperellipticCurve(f);

P := PeriodMatrix(X);
EndoRep := GeometricEndomorphismRepresentation(P, F);

print "";
print "Geometric endomorphism representations:";
print EndoRep;

idems := IsotypicalIdempotents(P, EndoRep);
idem := idems[1];

Q, proj := ComponentFromIdempotent(P, idem);
A, R := Explode(proj);

print "Induced action on homology:";
print R;

EQ := InducedPolarization(StandardSymplecticMatrix(3), R);
test, EQ := IsPolarization(EQ, Q);

Q, proj := ComponentFromIdempotent(P, idem : ProjOrInc := "Inc");
A, R := Explode(proj);

EQ := InducedPolarization(StandardSymplecticMatrix(3), R : ProjOrInc := "Inc");
test, EQ := IsPolarization(EQ, Q);

print "Sanity check for induced polarization:";
print test;
print EQ;
print PolarizationBasis(Q);

EQ0, T := FrobeniusFormAlternatingAlt(EQ);
print "Check claim in documentation:";
print EQ0 eq T*EQ*Transpose(T);

print "EQ0:";
print EQ0;
print BaseRing(EQ0);

Us := IsogenousPPLattices(EQ);
print "Isogenous lattices:";
print Us;

print "Corresponding determinants and polarizations:";
for U in Us do
    print Determinant(U);
    print U*EQ*Transpose(U);
end for;

print "Checking isogenous lattices via cover...";
Us := IsogenousPPLattices(EQ);
for U in Us do
    Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
    Qnew := Q*ChangeRing(Transpose(U), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
end for;
print "done";
