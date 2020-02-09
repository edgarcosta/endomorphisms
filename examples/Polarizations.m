SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

z12 := Exp(2*CC.1*Pi(CC)/12);
P := Matrix(CC, [
[z12^0,z12,z12^2,z12^3,z12^4,z12^5,z12^6,z12^7],
[z12^0,z12^2,z12^4,z12^6,z12^8,z12^10,z12^12,z12^14],
[z12^0,z12^4,z12^8,z12^12,z12^16,z12^20,z12^24,z12^28],
[z12^0,z12^7,z12^14,z12^21,z12^28,z12^35,z12^42,z12^49]
]);

test, E := SomePrincipalPolarization(P);
_, U := FrobeniusFormAlternatingRational(E);
Q := P*ChangeRing(Transpose(U), BaseRing(P));

/* This one also works over a usual ComplexField */
print SymplecticAutomorphismsCC(Q);
print SymplecticAutomorphisms(Q, F);
