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


intrinsic RosatiInvolution(EndoRep::SeqEnum, AorR::AlgMatElt) -> AlgMatElt
{Returns the Rosati involution of a tangent or homology representation. Assumes that the source has the standard principal polarization.}

As := [ gen[1] : gen in EndoRep ]; Rs := [ gen[2] : gen in EndoRep ];
g := #Rows(As[1]); IsR := #Rows(AorR) eq 2*g;

if IsR then
    R := AorR;
else;
    A := AorR;
    if not IsExact(Parent(A)) then
        Error("Exact input needed"); return 0;
        //s := Eltseq(MatrixInBasis(AorR, gensApp));
        //s := [ Round(c) : c in s ];
    end if;
    test, s := MatrixInBasis(A, As);
    assert test;
    s := Eltseq(s);
    R := &+[ s[i] * Rs[i] : i in [1..#Rs] ];
end if;

J := StandardSymplecticMatrix(g);
Rdagger := -J * Transpose(R) * J;
if IsR then
    return Rdagger;
else
    test, sdagger := MatrixInBasis(Rdagger, Rs);
    assert test;
    sdagger := Eltseq(sdagger);
    Adagger := &+[ sdagger[i] * As[i] : i in [1..#Rs] ];
    return Adagger;
end if;

end intrinsic;


intrinsic DegreeEstimate(EndoRep::SeqEnum, A::AlgMatElt) -> RngIntElt
{Estimates degree of corresponding endomorphism.}

Adagger := RosatiInvolution(EndoRep, A);
tr := Trace(A * Adagger) * Factorial(#Rows(A) - 1);
if IsExact(Parent(A)) then
    return (Integers() ! tr);
else
    Error("Exact input needed"); return 0;
end if;

end intrinsic;


intrinsic RosatiFixedModule(EndoRep::SeqEnum) -> SeqEnum
{Gives a basis of the ZZ-module of homological representations that are fixed
under Rosati.}

As := [ gen[1] : gen in EndoRep ]; Rs := [ gen[2] : gen in EndoRep ];
J := StandardSymplecticMatrix(#Rows(Rs[1]) div 2);
Rdiffs := [ (-J * Transpose(R) * J) - R : R in Rs ];
rowsM := [ ];
for Rdiff in Rdiffs do
    test, s := MatrixInBasis(Rdiff, Rs);
    assert test;
    Append(~rowsM, Eltseq(s));
end for;
M := Matrix(rowsM);
B := Basis(Kernel(M));
Asfixed := [ &+[ b[i]*As[i] : i in [1..#Rs] ] : b in B ];
Rsfixed := [ &+[ b[i]*Rs[i] : i in [1..#Rs] ] : b in B ];
EndoRepFixed := [ [* Asfixed[i], Rsfixed[i] *] : i in [1..#Rsfixed] ];
return EndoRepFixed;

end intrinsic;
