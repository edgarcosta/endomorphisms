/***
 *  Automorphisms
 *
 *  Copyright (C) 2016-2019
 *            Nils Bruin       (nbruin@sfu.ca)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
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
