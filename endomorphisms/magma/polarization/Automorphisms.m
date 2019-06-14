/***
 *  Automorphisms
 *
 *  Copyright (C) 2016-2019
 *            Nils Bruin       (nbruin@sfu.ca)
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
    if -E*Transpose(R)*E*R eq IdentityMatrix(Integers(), 2*g) then
        Append(~auts, R); Append(~auts, -R);
    end if;
end for;

G := MatrixGroup< 2*g, Integers() | auts >;
min := G ! (-G ! 1);
Q := quo< G | min >;
return Q, auts;

end intrinsic;
