/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Hanselman  (jeroen.hanselman@uni-ulm.de)
 *            Jeroen Sijsling   (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function MatrixExtra(P)
// Creates matrix over relevant base ring
CC := BaseRing(P);
if not assigned CC`epscomp then
    CC := ComplexFieldExtra(Precision(BaseRing(P)));
end if;
return ChangeRing(P, CC);

end function;


intrinsic StandardSymplecticMatrix(g::RngIntElt) -> AlgMatElt
{Standard symplectic 2 g x 2 g matrix.}

A := ScalarMatrix(g, 0); B := ScalarMatrix(g, 1); C := -B; D := A;
E := VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));
return ChangeRing(E, Rationals());

end intrinsic;

intrinsic AdjugateMatrix(M::AlgMatElt) -> AlgMatElt
{The adjugate matrix, i.e., Adj(A)*A = det(A) = A*Adj(A)}
    return Matrix(BaseRing(M), [[Cofactor(M, j, i) : j in [1..Ncols(M)]] : i in [1..Nrows(M)]]);
end intrinsic;

intrinsic PrincipalMinors(M::ModMatFldElt) -> SeqEnum[FldElt]
{Return the principal minors of M}
    require Nrows(M) eq Ncols(M): "the matrix needs to be square";
    return [Minor(M, I, I) where I := [j..n] : j  in [1..n]] where n:=Nrows(M);
end intrinsic;


intrinsic PolarizationBasis(P::ModMatFldElt) -> SeqEnum
{Determines a basis of the alternating forms giving rise to a polarization on the period matrix P.}

P := MatrixExtra(P);
JP := ComplexStructure(P); RR := BaseRing(JP);
gP := #Rows(JP) div 2; n := 4*gP^2;

/* Building a formal matrix corresponding to all possible polarisations */
R := PolynomialRing(RR, n); vars := GeneratorsSequence(R);
M := Matrix(R, 2*gP, 2*gP, vars);
JP_R := ChangeRing(JP, R);

/* Conditions that ensure that E(ix,iy) = E(x,y) and that E is antisymmetric */
Comm := Eltseq(Transpose(JP_R) * M * JP_R - M) cat Eltseq(M + Transpose(M));

/* Splitting previous linear equations by formal variable */
M :=  Matrix(RR, [ [ MonomialCoefficient(c, var) : c in Comm ] : var in vars ]);
vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Calculating polarization basis...";
Ker := IntegralLeftKernel(M);

/* Culling the correct polarizations using the conditions on E */
RR := BaseRing(JP); Es := [];
for r in Rows(Ker) do
    ht := Max([ Height(c) : c in Eltseq(r) ]);
    E := Matrix(Rationals(), 2*gP, 2*gP, Eltseq(r));
    ERR := ChangeRing(E, RR);
    /* Culling the correct polarizations using the conditions on E */
    Comm1 := Transpose(JP) * ERR * JP - ERR;
    if &and([Abs(c) lt ht*RR`epscomp : c in Eltseq(Comm1)]) then
        Comm2 := ERR + Transpose(ERR);
        if &and([Abs(c) lt ht*RR`epscomp : c in Eltseq(Comm2)]) then
            Append(~Es, E);
        end if;
    end if;
end for;
vprint EndoFind, 2 : "done calculating polarization basis.";
return Es;

end intrinsic;


intrinsic IsPolarization(E::., P::.) -> .
{Tests whether or not E defines a polarization on P.}

P := MatrixExtra(P);
CC := BaseRing(P);
/* Test Hermitian property */
ECC := ChangeRing(E, CC);
prod1 := P * ECC^(-1) * Transpose(P);
if not Maximum([ Abs(c) : c in Eltseq(prod1) ]) lt 10^20*BaseRing(P)`epscomp then
    return false;
end if;

/* Test positivity */
TPc := Matrix([ [ Conjugate(c) : c in Eltseq(r) ] : r in Rows(Transpose(P)) ]);
prod2 := CC.1 * P * ECC^(-1) * TPc;
res := [ Re(ev[1]) : ev in Eigenvalues(prod2) ];
if &and([ re gt CC`epsinv : re in res ]) then
    return true;
else
    return false;
end if;

end intrinsic;


function RealPolynomial(f)
    coeff, mon := CoefficientsAndMonomials(f);
    return &+[Real(c)*mon[i] : i->c in coeff];
end function;



function PfaffianAndMinors(Es, P)
    CC<I> := BaseRing(P);
    n := #Es;
    R<[x]> := PolynomialRing(Integers(), n);
    RCC<[y]> := PolynomialRing(CC, n);
    ER := &+[R.i * ChangeRing(Ei, R) : i->Ei in Es];
    pf := Pfaffian(ER);
    PRCC := ChangeRing(P, RCC);
    PRCC_star := ChangeRing(Conjugate(Transpose(P)), RCC);
    ERCC_inverse := ChangeRing(AdjugateMatrix(ER), RCC);
    M := I*PRCC*ERCC_inverse*PRCC_star;
    minors := PrincipalMinors(M);
    return pf, [RealPolynomial(elt) : elt in minors];
end function;

function _IsPolarization(coord, pf, minors : Principal:=false)
// Tests whether or not &+[ coord[i]*Es[i] : i in [1..n] ] defines a polarization on P.

// check determinant
det := Evaluate(pf, coord)^2;
if det eq 0 or (Principal and det ne 1) then
    return false;
end if;
// check positivity
CC := BaseRing(Parent(minors[1]));
for m in minors do
    if Evaluate(m, coord) lt CC`epsinv then
        return false;
    end if;
end for;

return true;
end function;


intrinsic SomePolarization(P::ModMatFldElt : B:=2, Principal:=false) -> AlgMatElt
{Tries to return some polarization for P.}

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P); n := #Es;
pf, minors := PfaffianAndMinors(Es, P);

counter := 0;
while true do
    counter +:= 1;
    test_power, exp := IsPowerOf(counter, 2);
    if test_power and exp ge 15 then
        B +:= 1;
    end if;
    if test_power and exp eq 25 then
        return false, 0;
    end if;
    D := [ -B..B ];

    CP := CartesianPower(D, n);
    tup := Random(CP);
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if _IsPolarization(tup, pf, minors : Principal:=Principal) then
        return true, &+[ tup[i]*Es[i] : i in [1..n] ];
    end if;
end while;

end intrinsic;

intrinsic SomePrincipalPolarization(P::ModMatFldElt : B := 2) -> AlgMatElt
{Tries to return some principal polarization for P.}
    return SomePolarization(P : B:=B, Principal:=true);
end intrinsic;

intrinsic SomePolarizations(P::ModMatFldElt : D:=[-2..2], Principal:=false) -> SeqEnum
{Tries to return some principal polarization for P.}
/* TODO: Implement Narasimhan--Nori if possible */

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P);
n := #Es; CP := CartesianPower(D, n);
pf, minors := PfaffianAndMinors(Es, P);


Es0 := [ ];
for tup in CP do
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if _IsPolarization(tup, pf, minors : Principal:=Principal) then
        if E in Es0 then
            Append(~Es0, E);
        end if;
    end if;
end for;
return Es0;

end intrinsic;

intrinsic SomePrincipalPolarizations(P::ModMatFldElt : B := 2) -> AlgMatElt
{Tries to return some principal polarization for P.}
    return SomePolarizations(P : B:=B, Principal:=true);
end intrinsic;


intrinsic SomePolarizationOld(P::ModMatFldElt : B := 2) -> AlgMatElt
{Tries to return some polarization for P.}

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P); n := #Es;

counter := 0;
while true do
    counter +:= 1;
    test_power, exp := IsPowerOf(counter, 2);
    if test_power and exp ge 15 then
        B +:= 1;
    end if;
    if test_power and exp eq 25 then
        return false, 0;
    end if;
    D := [ -B..B ];

    CP := CartesianPower(D, n);
    tup := Random(CP);
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Determinant(E) ne 0 then
        test := IsPolarization(E, P);
        if test then
            return true, E;
        end if;
    end if;
end while;
end intrinsic;

intrinsic SomePrincipalPolarizationOld(P::ModMatFldElt : B := 2) -> AlgMatElt
{Tries to return some principal polarization for P.}

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P); n := #Es;

counter := 0;
while true do
    counter +:= 1;
    test_power, exp := IsPowerOf(counter, 2);
    if test_power and exp ge 15 and IsEven(exp) then
        B +:= 1;
    end if;
    if test_power and exp eq 25 then
        return false, 0;
    end if;
    D := [ -B..B ];

    CP := CartesianPower(D, n);
    tup := Random(CP);
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Abs(Determinant(E)) eq 1 then
        test := IsPolarization(E, P);
        if test then
            return true, E;
        end if;
    end if;
end while;

end intrinsic;



intrinsic SomePolarizationsOld(P::ModMatFldElt : D := [-2..2]) -> SeqEnum
{Tries to return some principal polarization for P.}
/* TODO: Implement Narasimhan--Nori if possible */

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P);
n := #Es; CP := CartesianPower(D, n);

Es0 := [ ];
for tup in CP do
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Determinant(E) ne 0 then
        test := IsPolarization(E, P);
        if test and not E in Es0 then
            Append(~Es0, E);
        end if;
    end if;
end for;
return Es0;

end intrinsic;


intrinsic SomePrincipalPolarizationsOld(P::ModMatFldElt : D := [-2..2]) -> SeqEnum
{Tries to return some principal polarization for P.}
/* TODO: Implement Narasimhan--Nori if possible */

P := MatrixExtra(P);
Es := PolarizationBasis(P); CC := BaseRing(P);
n := #Es; CP := CartesianPower(D, n);

Es0 := [ ];
for tup in CP do
    E := &+[ tup[i]*Es[i] : i in [1..n] ];
    if Abs(Determinant(E)) eq 1 then
        test := IsPolarization(E, P);
        if test and not E in Es0 then
            Append(~Es0, E);
        end if;
    end if;
end for;
return Es0;

end intrinsic;
