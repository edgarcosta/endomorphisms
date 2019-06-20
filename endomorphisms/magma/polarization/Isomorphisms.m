/***
 *  Automorphisms
 *
 *  Copyright (C) 2016-2019
 *            Nils Bruin       (nbruin@sfu.ca)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

intrinsic SymplecticIsomorphismsCC(P::ModMatFldElt, Q::ModMatFldElt : GeoHomRep := [ ]) -> .
{Determines symplectic isomorphisms between abelian varieties defined by P and Q, which are assumed to be equipped with the standard principal polarization.}

if #Rows(P) ne #Rows(Q) then
    return [ ];
end if;

if #GeoHomRep eq 0 then
    GeoHomRep := GeometricHomomorphismRepresentationCC(P, Q);
end if;
Rs := [ tup[2] : tup in GeoHomRep ];

GeoHomRep := GeometricHomomorphismRepresentationCC(P, Q);

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

auts := SymplecticIsomorphismsCC(P, P : GeoHomRep := GeoHomRep);
Rs := [ aut[2] : aut in auts ];
g := #Rows(P);
GM := MatrixGroup< 2*g, Integers() | Rs >;
min := GM ! (-GM ! 1);
Q := quo< GM | min >;
G := Image(PermutationRepresentation(GM));
return G, Q, GM;

end intrinsic;
