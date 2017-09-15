/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Hanselman  (jeroen.hanselman@uni-ulm.de)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic StandardSymplecticMatrix(g::RngIntElt) -> .
{Standard symplectic 2 g x 2 g matrix.}

A := ScalarMatrix(g, 0); B := ScalarMatrix(g, 1); C := -B; D := A;
return VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));

end intrinsic;


function IntegralSymplecticBasisRecursive(E)
/* {Given an alternating integral matrix E, outputs an integral matrix T such
 * that T * E * Transpose(T) is in alternative symplectic standard form.} */

/* Use the Smith normal form to find two elements with smallest possible
 * pairing value and that can be extended to a basis: */
d := #Rows(E);
S, P, Q := SmithForm(E);
v1 := Matrix([ Eltseq(Rows(P)[1]) ]);
v2 := Matrix([ Eltseq(Rows(Transpose(Q))[1]) ]);
top := VerticalJoin(v1, v2);
if d eq 2 then
    return top;
end if;
K := Kernel(Transpose(VerticalJoin(v1*E, v2*E)));
B := Basis(K);
F := Matrix(Integers(), [ [ (Matrix(B[i]) * E * Transpose(Matrix(B[j])))[1,1] : j in [1..#B] ] : i in [1..#B] ]);
FSym := IntegralSymplecticBasisRecursive(F);
bottom := FSym * Matrix([ b : b in B ]);
T := VerticalJoin(top, bottom);
return T;

end function;


intrinsic IntegralSymplecticBasis(E::AlgMatElt) -> AlgMatElt
{Given an alternating integral matrix E, outputs an integral matrix T such that T * E * Transpose(T) is in symplectic standard form.}

d := #Rows(E);
rows := Rows(IntegralSymplecticBasisRecursive(E));
indices := [ i : i in [1..d] | IsOdd(i) ] cat [ i : i in [1..d] | IsEven(i) ];
T := Matrix(Integers(), [ [ c : c in Eltseq(rows[index]) ] : index in indices ]);
return T;

end intrinsic;


intrinsic PrincipalBasis(E::AlgMatElt) -> AlgMatElt
{Given an alternating integral matrix E, outputs a matrix T such that T * E * Transpose(T) is in principal standard form.}

T := IntegralSymplecticBasis(E);
T := ChangeRing(T, Rationals());
S := T * E * Transpose(T);
g := #Rows(E) div 2;
rows := Rows(T);
rows1 := rows[1..g];
rows2 := [ (1/S[i, g + i])*rows[g + i] : i in [1..g] ];
T := Matrix(Rationals(), [ [ c : c in Eltseq(row) ] : row in rows1 cat rows2 ]);
return T;

end intrinsic;


intrinsic FindPolarizationBasis(P::.) -> .
{Determines a basis of the alternating forms giving rise to a polarization on the period matrix P.}

JP :=ComplexStructure(P); RR := BaseRing(JP);
gP := #Rows(JP) div 2; n := 4 * gP ^2;

/* Building a formal matrix corresponding to all possible polarisations */
R := PolynomialRing(RR, n); vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gP, 2 *gP, vars);
JP_R := ChangeRing(JP, R);

/* Conditions that ensure that E(ix,iy) = E(x,y) and that E is anti-symmetric */
Comm := Eltseq(JP_R* M * Transpose(JP_R)- M) cat Eltseq(M+Transpose(M));

/* Splitting previous linear equations by formal variable */
M :=  Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] :var in vars ]);
Ker := IntegralLeftKernel(M);

/* Culling the correct polarizations using the conditions on E */
RR:=BaseRing(JP); Pols := [];
for r in Rows(Ker) do
    alpha := Matrix(Rationals(), 2*gP, 2*gP, Eltseq(r));
    /* Culling the correct polarizations using the conditions on E */
    Comm := JP*alpha*Transpose(JP) - alpha;
    Comm2 := alpha + Transpose(alpha);
    if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm)]) then
        if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm2)]) then
            Append(~Pols, alpha);
        end if;
    end if;
end for;
return Pols;

end intrinsic;


/* TODO: This function needs integrality properties */
intrinsic FindSymplecticBasis(M::.) -> .
{Determines a symplectic basis of the module M.}

V := SymplecticSpace(Matrix(M));
S := HyperbolicSplitting(V);
B := &cat[ [ Vector(v) : v in S[1][i] ] : i in [1..n] ];
return Matrix(B);
end intrinsic;


/* TODO: This function needs integrality properties and should return a
 * principal polarization on P itself if one exists */
intrinsic FindPrincipalPolarization(P::.) -> .
{Finds a principal por}

M := FindPolarizationBasis(P)[1];
N := FindSymplecticBasis(M);
N *:= &*[ Denominator(c) : c in Eltseq(N) ];
return N;

end intrinsic;
