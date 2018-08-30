AttachSpec("../../endomorphisms/magma/spec");

prec := 100;
R<x> := PolynomialRing(Rationals());
CC := ComplexFieldExtra(prec);

f := x*(x^4 + 2*x^3 - 3*x^2 + 7*x - 1);
g := x^8 + 2*x^6 - 3*x^4 + 7*x^2 - 1;
X := HyperellipticCurve(f);
Y := HyperellipticCurve(g);

eqsXCC := EmbedCurveEquations(X, prec);
eqsXF := DefiningEquations(X);
P := PeriodMatrix(eqsXCC, eqsXF);
eqsYCC := EmbedCurveEquations(Y, prec);
eqsYF := DefiningEquations(Y);
Q := PeriodMatrix(eqsYCC, eqsYF);

P := ChangeRing(P, CC); Q := ChangeRing(Q, CC);
GeoHomoRep := GeometricHomomorphismRepresentationPartial(P, Q);
print "Number of elements in isogeny basis:", #GeoHomoRep;
A, R := Explode(GeoHomoRep[1]);
RCC := ChangeRing(R, CC);
comm := A*P - Q*RCC;

print "Test almost 0:";
print Maximum([ Abs(c) : c in Eltseq(comm) ]);

print "Rank of R:";
print Rank(R);

/* TODO: Change these algorithms (for Jeroen) */

h := [* A, R *];
print "Connected component group of kernel:";
KCCG := KerModKer0(h, P, Q);
print KCCG;

print "Trivial component of kernel:";
K, inc := Ker0(h, P, Q);
itan, ihom := Explode(inc);
if K ne 0 then
    print ChangeRing(K, ComplexField(5));
    print "Inclusion on tangent:";
    print ChangeRing(itan, ComplexField(5));
    print "Inclusion on homology:";
    print ihom;
    print "Test almost 0:";
    ihomCC := ChangeRing(ihom, CC);
    comm := K*itan - ihomCC*P;
    print Maximum([ Abs(c) : c in Eltseq(comm) ]);
else
    print "trivial";
end if;

print "Cokernel:";
C, proj := Coker(h, P, Q);
ptan, phom := Explode(proj);
if C ne 0 then
    print ChangeRing(C, ComplexField(5));
    print "Projection on tangent:";
    print ChangeRing(ptan, ComplexField(5));
    print "Projection on homology:";
    print phom;
    print "Test almost 0:";
    phomCC := ChangeRing(phom, CC);
    comm := Q*ptan - phomCC*C;
    print Maximum([ Abs(c) : c in Eltseq(comm) ]);
else
    print "trivial";
end if;

print "Image:";
I, inc := Img(h, P, Q);
itan, ihom := Explode(inc);
if I ne 0 then
    print ChangeRing(I, ComplexField(5));
    print "Inclusion on tangent:";
    print ChangeRing(itan, ComplexField(5));
    print "Inclusion on homology:";
    print ihom;
    print "Test almost 0:";
    ihomCC := ChangeRing(ihom, CC);
    comm := I*itan - ihomCC*Q;
    print Maximum([ Abs(c) : c in Eltseq(comm) ]);
else
    print "trivial";
end if;

exit;
