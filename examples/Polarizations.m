SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

/*
prec := 200;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;
*/
prec := 200;
CC := ComplexField(prec);

z12 := Exp(2*CC.1*Pi(CC)/12);
P := Matrix(CC, [
[z12^0,z12,z12^2,z12^3,z12^4,z12^5,z12^6,z12^7],
[z12^0,z12^2,z12^4,z12^6,z12^8,z12^10,z12^12,z12^14],
[z12^0,z12^4,z12^8,z12^12,z12^16,z12^20,z12^24,z12^28],
[z12^0,z12^7,z12^14,z12^21,z12^28,z12^35,z12^42,z12^49]
]);

Es := PolarizationBasis(P);
print Es;
print SomePrincipalPolarization(P);
