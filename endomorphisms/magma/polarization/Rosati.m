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


intrinsic RosatiInvolution(GeoEndoRep::SeqEnum, AorR::AlgMatElt) -> AlgMatElt
{Returns the Rosati involution of a tangent or homology representation.}

As := [ gen[1] : gen in GeoEndoRep ]; Rs := [ gen[2] : gen in GeoEndoRep ];
g := #Rows(As[1]); isR := #Rows(AorR) eq 2*g;

if IsR then
    R := AorR;
else;
    A := AorR;
    if not IsExact(Parent(A)) then
        Error("Exact input needed"); return 0;
        //s := Eltseq(MatrixInBasis(AorR, gensApp));
        //s := [ Round(c) : c in s ];
    end if;
    s := Eltseq(MatrixInBasis(A, As));
    R := &+[ s[i] * Rs[i] : i in [1..#Rs] ];
end if;

J := StandardSymplecticMatrix(g);
Rdagger := -J * Transpose(R) * J;
if IsR then
    return Rdagger;
else
    sdagger := Eltseq(MatrixInBasis(Rdagger, Rs));
    Adagger := &+[ sdagger[i] * B[i] : i in [1..#Rs] ];
    return Adagger;
end if;

end intrinsic;


intrinsic DegreeEstimate(GeoEndoRep::SeqEnum, A::AlgMatElt) -> RngIntElt
{Estimates degree of corresponding endomorphism.}

Adagger := RosatiInvolution(GeoEndoRep, A);
tr := Trace(A * Adagger) * Factorial(#Rows(A) - 1);
if IsExact(Parent(A)) then
    return (Integers() ! tr);
else
    Error("Exact input needed"); return 0;
    //return Round(tr);
end if;

end intrinsic;


intrinsic RosatiFixedModule(GeoEndoRep::SeqEnum) -> SeqEnum
{Gives a basis of the ZZ-module of homological representations that are fixed
under Rosati.}

Rs := [ gen[2] : gen in GeoEndoRep ];
J := StandardSymplecticMatrix(#Rows(Rs[1]) div 2);
Rdiffs := [ (-J * Transpose(R) * J) - R : R in Rs ];
M := Matrix([ Eltseq(MatrixInBasis(Rdiff, Rs)) : Rdiff in Rdiffs ]);
B := Basis(Kernel(M));
Rsfixed := [ &+[ b[i]*Rs[i] : i in [1..#Rs] ] : b in B ];
return Rsfixed;

end intrinsic;
