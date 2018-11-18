/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Hanselman  (jeroen.hanselman@uni-ulm.de)
 *            Jeroen Sijsling   (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic StandardSymplecticMatrix(g::RngIntElt) -> AlgMatElt
{Standard symplectic 2 g x 2 g matrix.}

A := ScalarMatrix(g, 0); B := ScalarMatrix(g, 1); C := -B; D := A;
E := VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));
return ChangeRing(E, Rationals());

end intrinsic;


intrinsic PolarizationBasis(P::ModMatFldElt) -> SeqEnum
{Determines a basis of the alternating forms giving rise to a polarization on the period matrix P.}

JP := ComplexStructure(P); RR := BaseRing(JP);
gP := #Rows(JP) div 2; n := 4 * gP^2;

/* Building a formal matrix corresponding to all possible polarisations */
R := PolynomialRing(RR, n); vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gP, 2 * gP, vars);
JP_R := ChangeRing(JP, R);

/* Conditions that ensure that E(ix,iy) = E(x,y) and that E is antisymmetric */
Comm := Eltseq(JP_R * M * Transpose(JP_R) - M) cat Eltseq(M + Transpose(M));

/* Splitting previous linear equations by formal variable */
M :=  Matrix(RR, [ [ MonomialCoefficient(c, var) : c in Comm ] : var in vars ]);
Ker := IntegralLeftKernel(M);

/* Culling the correct polarizations using the conditions on E */
RR := BaseRing(JP); Es := [];
for r in Rows(Ker) do
    E := Matrix(Rationals(), 2*gP, 2*gP, Eltseq(r));
    ERR := ChangeRing(E, RR);
    /* Culling the correct polarizations using the conditions on E */
    Comm1 := JP * ERR * Transpose(JP) - ERR;
    if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm1)]) then
        Comm2 := ERR + Transpose(ERR);
        if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm2)]) then
            Append(~Es, E);
        end if;
    end if;
end for;
return Es;

end intrinsic;


intrinsic IsPolarization(E::., P::.) -> .
{Tests whether or not E defines a polarization on P.}

CC := BaseRing(P);
/* Test Hermitian property */
ECC := ChangeRing(E, CC);
prod1 := P * ECC * Transpose(P);
if not Maximum([ Abs(c) : c in Eltseq(prod1) ]) lt BaseRing(P)`epscomp then
    return false, 0;
end if;

/* Test positivity */
TPc := Matrix([ [ Conjugate(c) : c in Eltseq(r) ] : r in Rows(Transpose(P)) ]);
prod2 := BaseRing(P).1 * P * ECC * TPc;
res := [ Re(ev[1]) : ev in Eigenvalues(prod2) ];
if &and([ re gt CC`epsinv : re in res ]) then
    return true, E;
elif &and([ -re gt CC`epsinv : re in res ]) then
    return true, -E;
else
    return false, 0;
end if;

end intrinsic;


intrinsic SomePrincipalPolarizations(P::ModMatFldElt : D := [-2..2]) -> SeqEnum
{Tries to return some principal polarization for P.}
/* TODO: Implement Narasimhan--Nori */

Es := PolarizationBasis(P); CC := BaseRing(P);
n := #Es; CP := CartesianPower(D, n);

Es0 := [ ];
for tup in CP do
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Abs(Determinant(E)) eq 1 then
        test, E := IsPolarization(E, P);
        if test and not E in Es0 then
            Append(~Es0, E);
        end if;
    end if;
end for;
return Es0;

end intrinsic;


/* Next two functions should be redundant */
function SinglePrincipallyPolarizedCover(P, E)
// Given a period matrix admitting a polarization by E, find a (in general
// non-trivial) cover of P with the property on which E induces a principal
// polarization.

E0, T := FrobeniusFormAlternating(ChangeRing(E, Integers()));
T := ChangeRing(T, Rationals());
/* Now T*E*Transpose(T) = E0 */
g := #Rows(E0) div 2;
rows := Rows(T);
rows1 := rows[1..g];
rows2 := [ (1/E0[i, g + i])*rows[g + i] : i in [1..g] ];
U := Matrix(Rationals(), [ [ c : c in Eltseq(row) ] : row in rows1 cat rows2 ]);
Ui := U^(-1);
Q := P*ChangeRing(Ui, BaseRing(P));
return Q, Ui;

end function;


function SinglePrincipallyPolarizedQuotient(P, E)
// Given a period matrix admitting a polarization by E, find a (in general
// non-trivial) quotient of P with the property on which E induces a principal
// polarization.

E0, T := FrobeniusFormAlternating(ChangeRing(E, Integers()));
T := ChangeRing(T, Rationals());
/* Now T*E*Transpose(T) = E0 */
g := #Rows(E0) div 2;
rows := Rows(T);
rows1 := rows[1..g];
rows2 := [ (1/E0[i, g + i])*rows[g + i] : i in [1..g] ];
lcm := LCM([ Integers() ! E0[i, g + i] : i in [1..g] ]);
U := lcm*Matrix(Rationals(), [ [ c : c in Eltseq(row) ] : row in rows1 cat rows2 ]);
Ui := U^(-1);
Q := P*ChangeRing(Ui, BaseRing(P));
return Q, Ui;

end function;
