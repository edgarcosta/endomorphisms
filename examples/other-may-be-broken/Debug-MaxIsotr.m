/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 0);

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

EndoAlg, EndoDesc := EndomorphismStructure(EndoRep);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];

/* Find period matrix of fixed */
idems := IdempotentsFromStructure(EndoData);
idem := idems[1];
Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);

EQ := InducedPolarization(StandardSymplecticMatrix(3), R);
EQ0, T := FrobeniusFormAlternatingAlt(EQ);

print "Check claim in documentation:";
print EQ0 eq T*EQ*Transpose(T);

v1 := Matrix(Rationals(), [[0,0,1,0]]);
v2 := Matrix(Rationals(), [[0,0,0,1]]);
EQ := ChangeRing(EQ, Rationals());
T := ChangeRing(T, Rationals());

print "EQ0:";
print EQ0;

print "Number of symplectic submodules at 2:";
print #SymplecticSubmodules(4, 2);
print #SymplecticSubmodules(2, 2);
print #SymplecticSubmodules(2, 4);
//print #SymplecticSubmodules(2, 6);

Us := IsogenousPPLatticesG2(EQ);
print "Isogenous lattices:";
print Us;

print "Resulting polarizations:";
for U in Us do
    print Determinant(U);
    print U*EQ*Transpose(U);
end for;

exit;


d := 6;
EQ := Matrix(QQ, [[0,0,d,0],[0,0,0,1],[-d,0,0,0],[0,-1,0,0]]);
Us := IsogenousPPLatticesG2(EQ);
for U in Us do
    if not U*EQ*Transpose(U) eq d*ChangeRing(StandardSymplecticMatrix(2), Rationals()) then
        print "Problem";
        print U*EQ*Transpose(U);
    end if;
end for;

