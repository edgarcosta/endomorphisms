/***
 *  Automorphisms
 *
 *  Copyright (C) 2016-2017
 *            Nils Bruin       (nbruin@sfu.ca)
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

intrinsic SymplecticIsomorphisms(P::., Q::. : GeoHomRep := [ ]) -> .
{Determines symplectic isomorphisms between abelian varieties defined by P and Q, which are assumed to be equipped with the standard polarization.}

if #Rows(P) ne #Rows(Q) then
    return [];
end if;
g := #Rows(P);

if #GeoHomRep eq 0 then
    GeoHomRep := GeometricHomomorphismRepresentationCC(P, Q);
end if;
Rs := [ tup[2] : tup in GeoHomRep ];
E := StandardSymplecticMatrix(g);

n := #Rs;
T := PolynomialRing(Integers(), n);
RsT := [ ChangeRing(R, T) : R in Rs ];
ET := ChangeRing(E, T);
RUniv := &+[ (T.i)*RsT[i] : i in [1..n] ];
trf := Trace(-ET*Transpose(RUniv)*ET*RUniv);
L := LatticeWithGram(SymmetricMatrix(trf));
tups := ShortVectors(L, 2*g : Proof := true);

auts := [ ];
for tup in tups do
    seq := Eltseq(tup[1]);
    R := &+[ seq[i]*Rs[i] : i in [1..n] ];
    Append(~auts, R); Append(~auts, -R);
end for;
G := MatrixGroup< 2*g, Integers() | auts >;
min := G ! (-G ! 1);
Q := quo< G | min >;
return auts, Q;

end intrinsic;
