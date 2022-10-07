/***
 *  Automorphisms
 *
 *  Copyright (C) 2016-2019
 *            Nils Bruin       (nbruin@sfu.ca)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *  Copyright (C) 2022
 *            Edgar Costa      (edgarc@mit.edu)
 *
 *  See LICENSE.txt for license details.
 */

import "Polarizations.m": MatrixExtra;


intrinsic SymplecticIsomorphismsCC(P::ModMatFldElt, Q::ModMatFldElt : GeoHomRep := [ ]) -> .
{Determines symplectic isomorphisms between abelian varieties defined by P and Q, which are assumed to be equipped with the standard principal polarization.}

P := MatrixExtra(P); Q := MatrixExtra(Q);
if #Rows(P) ne #Rows(Q) then
    return [ ];
end if;

if #GeoHomRep eq 0 then
    GeoHomRep := GeometricHomomorphismRepresentationCC(P, Q);
end if;
if #GeoHomRep eq 0 then
    return [ ];
end if;
Rs := [ tup[2] : tup in GeoHomRep ];

g := #Rows(P); CC := BaseRing(P);
E := StandardSymplecticMatrix(g);
As := [ gen[1] : gen in GeoHomRep ];
Rs := [ gen[2] : gen in GeoHomRep ];

S := PolynomialRing(Rationals(), #Rs);
R := &+[ S.i*Rs[i] : i in [1..#Rs] ];
Rdagger := -E * Transpose(R) * E;
tr := Trace(R*Rdagger);

M := Matrix(Rationals(), #Rs,#Rs, [ [ Derivative(Derivative(tr, i), j) : j in [1..#Rs] ] : i in [1..#Rs] ]);
Lat := LatticeWithGram(M);
tups := ShortVectors(Lat, 4*g : Proof := true);
vs := [ tup[1] : tup in tups ];
vs := vs cat [ -v : v in vs ];

isos := [ ];
for v in vs do
    seq := Eltseq(v);
    A := &+[ seq[i]*As[i] : i in [1..#Rs] ];
    R := &+[ seq[i]*Rs[i] : i in [1..#Rs] ];
    Rdagger := -E * Transpose(R) * E;
    if R*Rdagger eq 1 then
        Append(~isos, [* A, R *]);
    end if;
end for;
return isos;

end intrinsic;


intrinsic SymplecticAutomorphismsCC(P::ModMatFldElt : GeoHomRep := [ ]) -> .
{Determines symplectic automorphism group of the abelian variety defined by P, which is assumed to be equipped with the standard principal polarization.}

P := MatrixExtra(P);
auts := SymplecticIsomorphismsCC(P, P : GeoHomRep := GeoHomRep);
Rs := [ aut[2] : aut in auts ];
g := #Rows(P);
GM := MatrixGroup< 2*g, Integers() | Rs >;
min := GM ! (-GM ! 1);
G := Image(PermutationRepresentation(GM));
Q := quo< GM | min >;
return G, Q, GM;

end intrinsic;


intrinsic SymplecticIsomorphisms(P::ModMatFldElt, Q::ModMatFldElt, F::Fld : GeoHomRep := [ ]) -> .
{Determines symplectic isomorphisms between abelian varieties defined by P and Q, which are assumed to be equipped with the standard principal polarization.}

P := MatrixExtra(P); Q := MatrixExtra(Q);
isosPart := SymplecticIsomorphismsCC(P, Q : GeoHomRep := GeoHomRep);
isosPart := [ [* ChangeRing(iso[1], F`CC), iso[2] *] : iso in isosPart ];
/* Determine minimal polynomials needed */
seqPart := &cat[ Eltseq(iso[1]) : iso in isosPart ];

/* Use splitting field instead of number field since the resulting field is
 * normal */
vprint EndoFind : "";
vprint EndoFind : "Finding number field defined by isomorphisms...";
K, seq, hFK := NumberFieldExtra(seqPart, F);
vprint EndoFind : "done finding number field defined by isomorphisms:";
vprint EndoFind : K;

assert #seq eq #seqPart;
if #seq eq 0 then
    return [ ], hFK;
end if;

r := #Rows(isosPart[1][1]); c := #Rows(Transpose(isosPart[1][1]));
As := [ Matrix(K, r, c, seq[((k - 1)*r*c + 1)..(k*r*c)]) : k in [1..#isosPart] ];
isos := [ [* As[k], isosPart[k][2] *] : k in [1..#isosPart] ];

/* Final check for correctness */
for i in [1..#isos] do
    abs := Max([ Abs(c) : c in Eltseq(EmbedMatrixExtra(isos[i][1]) - isosPart[i][1]) ]);
    assert abs lt 10^20*BaseRing(P)`epscomp;
end for;
return isos, hFK;

end intrinsic;


intrinsic SymplecticAutomorphisms(P::ModMatFldElt, F::Fld : GeoHomRep := [ ]) -> .
{Determines symplectic automorphism group of the abelian variety defined by P, which is assumed to be equipped with the standard principal polarization.}

P := MatrixExtra(P);
isosPart := SymplecticIsomorphismsCC(P, P : GeoHomRep := GeoHomRep);
isosPart := [ [* ChangeRing(iso[1], F`CC), iso[2] *] : iso in isosPart ];
/* Determine minimal polynomials needed */
seqPart := &cat[ Eltseq(iso[1]) : iso in isosPart ];

/* Use splitting field instead of number field since the resulting field is
 * normal */
vprint EndoFind : "";
vprint EndoFind : "Finding number field defined by isomorphisms...";
K, seq, hFK := SplittingFieldExtra(seqPart, F);
vprint EndoFind : "done finding number field defined by isomorphisms:";
vprint EndoFind : K;

assert #seq eq #seqPart;
if #seq eq 0 then
    return [ ], hFK;
end if;

r := #Rows(isosPart[1][1]); c := #Rows(Transpose(isosPart[1][1]));
As := [ Matrix(K, r, c, seq[((k - 1)*r*c + 1)..(k*r*c)]) : k in [1..#isosPart] ];
isos := [ [* As[k], isosPart[k][2] *] : k in [1..#isosPart] ];

/* Final check for correctness */
for i in [1..#isos] do
    abs := Max([ Abs(c) : c in Eltseq(EmbedMatrixExtra(isos[i][1]) - isosPart[i][1]) ]);
    assert abs lt 10^20*BaseRing(P)`epscomp;
end for;
return isos, hFK;

end intrinsic;



intrinsic GeometricAutomorphismGroupViaPeriods(C::CrvHyp : prec:=100) -> GrpPerm, Map, Map, Map
{ The automorphism group G of C as a permutation group.
  The second value is a map from G to the group of automorphisms of C.
  The third value is the map of the base field to L, where the automorphism are defined.
  The fourth value gives the action C_L x G -> C_L
}
// this follows very much AutomorphismGroup in Geometry/CrvHyp/transform.m
// but replaces IsGL2Equivalent with Sympletic Automorphisms

//Auxiliar functions that are available in Magma src files, but not exposed
// from Geometry/CrvHyp/transform.m
function TransformPolynomial(f,n,mat)
    // Transforms f, considered as a homogeneous polynomial of
    // degree n, by the matrix with entries mat = [a,b,c,d].
    a,b,c,d := Explode(mat);
    x1 := a*Parent(f).1 + b; z1 := c*Parent(f).1 + d;
    return &+[Coefficient(f,i)*x1^i*z1^(n-i) : i in [0..n]];
end function;

// from Geometry/CrvHyp/transform.m
function IsomorphismData(m)
    // Note that this function has the same name as that on elliptic
    // curves -- they should be merged into one function on MapIsoSch.

    C := Domain(m);
    P := CoordinateRing(Ambient(C));
    x := P.1; y := P.2; z := P.3;
    g := Genus(C);
    eqns := DefiningEquations(m);
    t := [ MonomialCoefficient(f,m) : m in [x,z], f in eqns[[1,3]] ];
    e := MonomialCoefficient(eqns[2], y);
    u := Polynomial(BaseRing(C), coeffs) where coeffs :=
        [ MonomialCoefficient(eqns[2],x^i*z^(g+1-i)) : i in [0..g+1] ];
    return t, e, u;
end function;

// from Geometry/CrvHyp/isomorphisms.m
function HyperellipticIsomorphism(C1,C2,t,e,u)
    P1 := CoordinateRing(Ambient(C1));
    x1 := P1.1; y1 := P1.2; z1 := P1.3;
    P2 := CoordinateRing(Ambient(C2));
    x2 := P2.1; y2 := P2.2; z2 := P2.3;
    g := Genus(C1);
    a,b,c,d := Explode(t);
    S := Eltseq(u);
    X2 := a*x1+b*z1; Z2 := c*x1+d*z1;
    Y2 := e*y1 + &+[ P1 | S[i+1]*x1^i*z1^(g+1-i) : i in [0..#S-1] ];
    det := a*d-b*c;
    X1 := d*x2-b*z2; Z1 := -c*x2+a*z2;
    Y1 := (1/e)*(det^(g+1)*y2
        - &+[ P2 | S[i+1]*X1^i*Z1^(g+1-i) : i in [0..#S-1] ]);
    /*
    print "[X2,Y2,Z2] =", [X2,Y2,Z2];
    print "[X1,Y1,Z1] =", [X1,Y1,Z1];
    */
return iso< C1 -> C2 | [ X2, Y2, Z2 ], [ X1, Y1, Z1 ] : Check:=false>;
end function;


g := Genus(C);
require g gt 0: "GeometricAutomorphismGroupViaPeriods not implemented for genus 0.";

K := BaseField(C);
// converts a curve to be over a NumberFieldExtra if it is not given in that way
if not assigned K`iota then
    C := CurveExtra(C : prec:=prec);
end if;
CC := Parent(K`iota);

C1, _ := SimplifiedModel(C);

vprintf CrvHypIso: "Computing period matrix...";
vtime CrvHypIso:
P := PeriodMatrix(C1);

vprintf CrvHypIso: "Computing symplectic automorphisms...";
vtime CrvHypIso:
list, hToL := SymplecticAutomorphisms(P, K);
vprintf CrvHypIso, 1: "SymplecticAutomorphisms:\n    %o\n", list;
GQbar := GeometricAutomorphismGroup(C1);
// check that we have the right number of elements
assert #list eq #GQbar;
L := Codomain(hToL);

// it is easier to create the map again than try to base change the map
CL := ChangeRing(C, L);
C1L, simp_map := SimplifiedModel(CL);
simp_map_inv := Inverse(simp_map);
f1 := HyperellipticPolynomials(C1L);
listA := [elt[1] : elt in list];
// move identity and involution to the front
for i->elt in [I, -I] where I is IdentityMatrix(L, 2) do
    j := Position(listA, elt);
    tmp := listA[i];
    listA[i] := elt;
    listA[j] := tmp;
end for;
// Construct the corresponding permutation group
G := PermutationGroup<#listA | [ [ Position(listA, elt1*elt0) : elt0 in listA ] : elt1 in listA ] >;
assert IdentifyGroup(G) eq IdentifyGroup(GQbar);

// Construct the isomorphisms
list1 := [ ];
for i->elt in listA do
    t := Reverse(Eltseq(elt));
    f2 := TransformPolynomial(f1, 6, t);
    c := LeadingCoefficient(f2)/LeadingCoefficient(f1);
    bool, rtc := IsSquare(c);
    assert bool;
    bool, iso := IsAutomorphism(HyperellipticIsomorphism(C1L, C1L, t, rtc, 0));
    assert bool;
    Append(~list1, simp_map * iso * simp_map_inv);
end for;

// The construction of the map is somewhat primitive, but should be OK,
// since the groups are small...
GT := car< G, Aut(CL) >;
m := map< G -> Aut(CL) |
          [ GT | <G.j, l> : j->l in list1 ] >;

return G, m, hToL, map< car<CL, G> -> CL | p :-> m(p[2])(p[1]) >;
end intrinsic;
