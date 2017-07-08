/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic StandardSymplecticMatrix(g::RngIntElt) -> .
{Standard symplectic 2 g x 2 g matrix.}

A := ScalarMatrix(g, 0); B := ScalarMatrix(g, -1); C := -B; D := A;
return VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));

end intrinsic;


intrinsic RosatiInvolution(GeoEndoRep::SeqEnum, A::.) -> .
{Returns the Rosati involution of A.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];
Rs := gensHom;
if IsExact(Parent(A)) then
    B := gensTan;
    s := Eltseq(MatrixInBasis(A, B));
else
    B := gensApp;
    s := Eltseq(MatrixInBasis(A, B));
    s := [ Round(c) : c in s ];
end if;
R := &+[ s[i] * Rs[i] : i in [1..#Rs] ];
J := StandardSymplecticMatrix(#Rows(A));
Rdagger := -J * Transpose(R) * J;
sdagger := Eltseq(MatrixInBasis(Rdagger, Rs));
Adagger := &+[ sdagger[i] * B[i] : i in [1..#Rs] ];
return Adagger;

end intrinsic;


intrinsic DegreeEstimate(GeoEndoRep::SeqEnum, A::.) -> .
{Estimates degree of corresponding endomorphism.}

Adagger := RosatiInvolution(GeoEndoRep, A);
tr := Trace(A * Adagger) * Factorial(#Rows(A) - 1);
if IsExact(Parent(A)) then
    return (Integers() ! tr);
else
    return Round(tr);
end if;

end intrinsic;


intrinsic RosatiFixedModule(GeoEndoRep::SeqEnum) -> .
{Gives the ZZ-module of homological representations that are fixed under Rosati.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ]; Rs := gensHom;
gensApp := [ gen[3] : gen in GeoEndoRep ];
J := StandardSymplecticMatrix(#Rows(Rs[1]) div 2);
Rdiffs := [ (-J * Transpose(rep[2]) * J) - rep[2] : rep in GeoEndoRep ];
M := Matrix([ Eltseq(MatrixInBasis(Rdiff, Rs)) : Rdiff in Rdiffs ]);
B := Basis(Kernel(M));
Rsfixed := [ &+[ b[i]*Rs[i] : i in [1..#Rs] ] : b in B ];
return Rsfixed;

end intrinsic;
