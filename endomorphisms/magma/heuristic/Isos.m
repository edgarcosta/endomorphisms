/***
 *  Finding isomorphisms
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic IsomorphismsCC(P::ModMatFldElt, Q::ModMatFldElt) -> SeqEnum
{Finds tangent representation over CC of isomorphisms.}

assert #Rows(P) eq #Rows(Q);
gens := GeometricHomomorphismRepresentationCC(P, Q);

g := #Rows(P); CC := BaseRing(P);
E := StandardSymplecticMatrix(g);
As := [ gen[1] : gen in gens ];
Rs := [ gen[2] : gen in gens ];

S := PolynomialRing(Rationals(), #Rs);
R := &+[ S.i*Rs[i] : i in [1..#Rs] ];
Rdagger := -E * Transpose(R) * E;
tr := Trace(R*Rdagger);

M := Matrix(Rationals(), #Rs,#Rs, [ [ Derivative(Derivative(tr, i), j) : j in [1..#Rs] ] : i in [1..#Rs] ]);
Lat := LatticeWithGram(M);
tups := ShortVectors(Lat, 4*g : Proof := true);
vs := [ tup[1] : tup in tups ];
vs := vs cat [ -v : v in vs ];

gens0 := [ ];
for v in vs do
    seq := Eltseq(v);
    A := &+[ seq[i]*As[i] : i in [1..#Rs] ];
    R := &+[ seq[i]*Rs[i] : i in [1..#Rs] ];
    Rdagger := -E * Transpose(R) * E;
    if R*Rdagger eq 1 then
        Append(~gens0, [* A, R *]);
    end if;
end for;

return gens0;

end intrinsic;
