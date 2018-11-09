//SetVerbose("EndoFind", 1);

prec := 600;
prec := 200;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := x^8 + x^6 + 2; h := R ! 0;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
Y := HyperellipticCurve(x^4 + x^3 + 2);
//Y := HyperellipticCurve(x*(x^4 + x^3 + 2));

P := PeriodMatrix(X);
Q := PeriodMatrix(Y);
HomRepPQ := GeometricHomomorphismRepresentationCC(P, Q);
HomRepQP := GeometricHomomorphismRepresentationCC(Q, P);
//HomRepPQ := GeometricHomomorphismRepresentation(P, Q, F);
//HomRepQP := GeometricHomomorphismRepresentation(Q, P, F);
hPQ := HomRepPQ[1];
hQP := HomRepQP[1];

print "Homomorphism:";
print hPQ[1];
print "Component group of kernel:";
print KerModKer0(hPQ, P, Q);
K, hK := Ker0(hPQ, P, Q);
print "Kernel:";
print hK[2];
C, hC := Coker(hPQ, P, Q);
print "Cokernel:";
print hC[2];
I, hI := ImgInc(hPQ, P, Q);
print "Image as inclusion:";
print hI[2];
I, hI := ImgProj(hPQ, P, Q);
print "Image as projection:";
print hI[2];

print "";
print "Homomorphism:";
print hQP[2];
print "Component group of kernel:";
print KerModKer0(hQP, Q, P);
K, hK := Ker0(hQP, Q, P);
print "Kernel:";
print hK[2];
C, hC := Coker(hQP, Q, P);
print "Cokernel:";
print hC[2];
I, hI := ImgInc(hQP, Q, P);
print "Image as inclusion:";
print hI[2];
I, hI := ImgProj(hQP, Q, P);
print "Image as projection:";
print hI[2];
