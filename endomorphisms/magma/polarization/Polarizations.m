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


/* TODO: Superseded by FrobeniusFormAlternating */
//function IntegralSymplecticBasisRecursive(E)
///* {Given an alternating integral matrix E, outputs an integral matrix T such
// * that T * E * Transpose(T) is in alternative symplectic standard form.} */
//
//E := Matrix(Integers(), [ [ c : c in Eltseq(row) ] : row in Rows(E) ]);
///* Use the Smith normal form to find two elements with smallest possible
// * pairing value and that can be extended to a basis: */
//d := #Rows(E);
//S, P, Q := SmithForm(E);
//v1 := Matrix([ Eltseq(Rows(P)[1]) ]);
//v2 := Matrix([ Eltseq(Rows(Transpose(Q))[1]) ]);
//top := VerticalJoin(v1, v2);
//if d eq 2 then
//    return top;
//end if;
//K := Kernel(Transpose(VerticalJoin(v1*E, v2*E)));
//B := Basis(K);
//F := Matrix(Integers(), [ [ (Matrix(B[i]) * E * Transpose(Matrix(B[j])))[1,1] : j in [1..#B] ] : i in [1..#B] ]);
//FSym := IntegralSymplecticBasisRecursive(F);
//bottom := FSym * Matrix([ b : b in B ]);
//T := VerticalJoin(top, bottom);
//return T;
//
//end function;


/* TODO: Superseded by FrobeniusFormAlternating */
//intrinsic IntegralSymplecticBasis(E::AlgMatElt) -> AlgMatElt
//{Given an alternating integral matrix E, outputs an integral matrix T such that T * E * Transpose(T) is in symplectic standard form.}
//
//d := #Rows(E);
//rows := Rows(IntegralSymplecticBasisRecursive(E));
//indices := [ i : i in [1..d] | IsOdd(i) ] cat [ i : i in [1..d] | IsEven(i) ];
//T := Matrix(Integers(), [ [ c : c in Eltseq(rows[index]) ] : index in indices ]);
//return T;
//
//end intrinsic;


/* TODO: This function needs integrality properties */
//intrinsic FindSymplecticBasis(M::.) -> .
//{Determines a symplectic basis of the module M.}
//
//V := SymplecticSpace(Matrix(M));
//S := HyperbolicSplitting(V);
//n := #S;
//B := &cat[ [ Vector(v) : v in S[1][i] ] : i in [1..n] ];
//return Matrix(B);
//end intrinsic;


intrinsic PrincipallyPolarizedCover(P::ModMatFldElt, E::AlgMatElt) -> ModMatFldElt
{Given a period matrix admitting a polarization by E, find a (in general non-trivial) quotient of P with the property on which E induces a principal polarization.}

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

end intrinsic;



intrinsic PrincipallyPolarizedQuotient(P::ModMatFldElt, E::AlgMatElt) -> ModMatFldElt
{Given a period matrix admitting a polarization by E, find a (in general non-trivial) quotient of P with the property on which E induces a principal polarization.}

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

end intrinsic;


intrinsic FindPolarizationBasis(P::ModMatFldElt) -> SeqEnum
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


intrinsic FindPrincipalPolarization(P::ModMatFldElt : D := [-2..2]) -> SeqEnum
{Finds some principal polarization for P.}
/* TODO: Implement Narasimhan--Nori */

Es := FindPolarizationBasis(P); CC := BaseRing(P);
n := #Es; CP := CartesianPower(D, n);

Es0 := [ ];
for tup in CP do
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Abs(Determinant(E)) eq 1 then
        /* Test Hermitian property */
        Ei := ChangeRing(E, CC);
        prod1 := P * Ei * Transpose(P);
        test1 := Maximum([ Abs(c) : c in Eltseq(prod1) ]) lt BaseRing(P)`epscomp;

        /* Test positivity */
        TPc := Matrix([ [ Conjugate(c) : c in Eltseq(r) ] : r in Rows(Transpose(P)) ]);
        prod2 := BaseRing(P).1 * P * Ei * TPc;
        res := [ Re(ev[1]) : ev in Eigenvalues(prod2) ];
        if &and([ re gt CC`epsinv : re in res ]) then
            test2 := true;
        elif &and([ -re gt CC`epsinv : re in res ]) then
            test2 := true;
            E *:= -1;
        else
            test2 := false;
        end if;
        if test1 and test2 then
            if not E in Es0 then
                Append(~Es0, E);
            end if;
        end if;
    end if;
end for;
return Es0;

end intrinsic;
