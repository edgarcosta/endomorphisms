SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 100;
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

z8 := Exp(2*CC.1*Pi(CC)/8);
P := Matrix(CC, [
[z8^0,z8^1,z8^2,z8^3,z8^4,z8^5],
[z8^0,z8^2,z8^4,z8^6,z8^8,z8^10],
[z8^0,z8^5,z8^10,z8^15,z8^20,z8^25]
]);

Q := Matrix(CC, [
[1,z12^4]
]);

Q := Matrix(CC, [
[1,z12^3]
]);

time GeoEndoRepCC := GeometricEndomorphismRepresentationCC(DiagonalJoin([ P, P ]));
time GeoEndoRepCC := GeometricEndomorphismRepresentationCC(DiagonalJoin([ P, P, P ]));
GeoEndoRepCC := GeometricEndomorphismRepresentationCC(P);
GeoHomRepCC := GeometricHomomorphismRepresentationCC(P, Q);
GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
print GeoEndoRep;
