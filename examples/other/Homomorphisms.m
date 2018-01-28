AttachSpec("../../endomorphisms/magma/spec");

CC := ComplexFieldExtra(1000); CCSmall := ComplexField(5);
R<x> := PolynomialRing(CC);
X := SE_Curve(x^6 + x^2 + 1, 2 : Prec := Precision(CC));
Y := SE_Curve(x^3 + x + 1, 2 : Prec := Precision(CC));
P := ChangeRing(X`BigPeriodMatrix, CC) / 2;
Q := ChangeRing(Y`BigPeriodMatrix, CC) / 2;

print "P to Q:";
GeoHomRep := GeometricHomomorphismRepresentationPartial(P, Q);
for tup in GeoHomRep do
    T := tup[1];
    R := tup[2];
    print "Check small:";
    print CCSmall ! Maximum([ Abs(c) : c in Eltseq(T*P - Q*ChangeRing(R, CC)) ]);
end for;
print "---";

Tn := TangentRepresentation(R, P, Q);
Rn := HomologyRepresentation(T, P, Q);

print "Check small:";
print CCSmall ! Maximum([ Abs(c) : c in Eltseq(Tn - T) ]);
print CCSmall ! Maximum([ Abs(c) : c in Eltseq(Rn - R) ]);
print "---";

print ChangeRing(T, CCSmall);
print R;
print "===";

print "Q to P:";
GeoHomRep := GeometricHomomorphismRepresentationPartial(Q, P);
for tup in GeoHomRep do
    T := tup[1];
    R := tup[2];
    print "Check small:";
    print CCSmall ! Maximum([ Abs(c) : c in Eltseq(T*Q - P*ChangeRing(R, CC)) ]);
end for;
print "---";

Tn := TangentRepresentation(R, Q, P);
Rn := HomologyRepresentation(T, Q, P);

print "Check small:";
print CCSmall ! Maximum([ Abs(c) : c in Eltseq(Tn - T) ]);
print CCSmall ! Maximum([ Abs(c) : c in Eltseq(Rn - R) ]);
print "---";

print ChangeRing(T, CCSmall);
print R;

exit;
