/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 0);

prec := 500;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x,y> := PolynomialRing(F, 2);
w := x*y + 12*y^2 + 32;
f := 108*w^2*y + 8208*w*y^2 + x^3 + 36*x^2*y + 96*x^2 - 6480*x*y^2 + 2304*x*y + 2928*x + 74304*y^3 - 248832*y^2 - 75456*y + 27584;
print f;

/*
w := x*y - 32*y^2 - 12;
f := 108*x^2 + 8208*x + w^3 + 96*w^2*y + 36*w^2 + 2928*w*y^2 + 2304*w*y - 6480*w + 27584*y^3 - 75456*y^2 - 248832*y + 74304;
print f;
*/

X := PlaneCurve(f);
print "Curve:";
print X;

P := PeriodMatrix(X);
print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
//L<s> := BaseRing(GeoEndoRep[1][1]);

print "";
print "Endomorphism representations:";
print GeoEndoRep;

lat, sthash := EndomorphismLattice(GeoEndoRep);
print "";
print "Endomorphism lattice:";
print lat;

rep := lat[2][1][2][1];
Rs := [ pair[2] : pair in rep ];
Rs0 := [ R : R in Rs | Determinant(R) eq 0 ];

R := Rs0[1]/2;
A := TangentRepresentation(R, P, P);
idem := [* A, R *];

/* This defines a valid numerical map */
Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);
print Max([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, BaseRing(P))) ]);

tau := Q[1,1]^(-1)*Q[1,2];
j := jInvariant(tau);

F := RationalsExtra(prec);
E := ReconstructCurveG1(Q, F);
print E;

K<r> := BaseField(E);
a4 := 1212597*r^4 + 989496*r;
a6 := -259785954*r^3 + 1111959306;
R<x> := PolynomialRing(K);
f := x^3 + a4*x + a6;
E := HyperellipticCurve(f);
Q := PeriodMatrix(E);

HomRepCC := GeometricHomomorphismRepresentationCC(P, Q);
HomRep := GeometricHomomorphismRepresentation(P, Q, K);

XK := ChangeRing(X, K);
// Weierstrass point:
P0 := XK ! [ -4*r^5 - 4*r^4 + 8*r^2 + 8*r - 32, 0, 1 ];

F := DefiningEquations(XK)[1];
S := Parent(F);
R := PolynomialRing(K);
B := 5; D := [-100..100];
CP := CartesianPower(D, 2);
for tup in [ tup : tup in CP | not (tup[1] eq 0) and (tup[2] eq 0) ] do
    h := hom< S -> R | [ R.1, tup[1], tup[2] ] >;
    h := hom< S -> R | [ tup[1], R.1, tup[2] ] >;
    h := hom< S -> R | [ tup[1], tup[2], R.1 ] >;
    if h(F) ne 0 then
        rts := [ rt[1] : rt in Roots(h(F)) ];
        if #rts ne 0 then
            for rt in rts do
                //P0 := XK ! [ rt, tup[1], tup[2] ];
                //P0 := XK ! [ tup[1], rt, tup[2] ];
                P0 := XK ! [ tup[1], tup[2], rt ];
                if not IsSingular(P0) then
                    if not IsWeierstrassPlace(Place(P0)) then
                        print P0;
                    end if;
                end if;
            end for;
        end if;
    end if;
end for;



/*
[ [*
    [1/91*(-38*r^3 + 22) 1/91*(-18*r^4 + 20*r) 1/91*(-216*r^4 - 1456*r^3 + 240*r
        + 728) 1/91*(-640*r^4 + 792*r)],

    [-1  0  1  0 -1  0  1  1]
    [ 1  0  0  0  0  0  1  1]
*], [*
    [1/91*(-22*r^3 - 16) 1/91*(2*r^4 + 18*r) 1/91*(24*r^4 - 728*r^3 + 216*r -
        728) 1/91*(152*r^4 + 640*r)],

    [-1 -1  1 -1 -1 -1  0  0]
    [ 0  1  1 -1  0  0  1  0]
*], [*
    [1/91*(-6*r^4 + 38*r^3 - 54*r - 22) 1/91*(18*r^4 - 38*r^3 - 20*r + 22)
        1/91*(-240*r^4 + 1000*r^3 - 2160*r - 464) 1/91*(640*r^4 - 1456*r^3 -
        792*r + 728)],

    [ 1 -1 -1 -1  0  0 -3 -3]
    [ 2 -1  0  2  2  0  0  0]
*], [*
    [1/91*(-54*r^4 + 60*r) 1/91*(16*r^3 - 38) 1/91*(-1920*r^4 + 192*r^3 + 2376*r
        - 456) 8*r^3 - 16],

    [ 0  0 -2  2  0 -1  0  1]
    [ 0  3  1 -1  0  2  0  1]
*] ]
*/

/*
Rs :=
[
Matrix([[-1,  0,  1,  0, -1,  0,  1,  1], [ 1,  0,  0,  0,  0,  0,  1,  1]]),
Matrix([[-1, -1,  1, -1, -1, -1,  0,  0], [ 0,  1,  1, -1,  0,  0,  1,  0]]),
Matrix([[ 1, -1, -1, -1,  0,  0, -3, -3], [ 2, -1,  0,  2,  2,  0,  0,  0]]),
Matrix([[ 0,  0, -2,  2,  0, -1,  0,  1], [ 0,  3,  1, -1,  0,  2,  0,  1]])
];
*/
