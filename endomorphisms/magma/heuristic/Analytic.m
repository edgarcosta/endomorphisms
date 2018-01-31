/***
 *  Finding an approximate basis for the geometric endomorphism ring
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic ComplexStructure(P::.) -> .
{Returns the complex structure that corresponds to the period matrix P. It is
the RR-linear transformation on homology that corresponds to multiplication
with i.}

CC := BaseRing(P); RR := RealField(CC);
PSplit := VerticalSplitMatrix(P); iPSplit := VerticalSplitMatrix(CC.1 * P);
return NumericalRightSolve(PSplit, iPSplit);

end intrinsic;


intrinsic RationalHomomorphismEquations(JP::., JQ::.) -> .
{Given two complex structures JP and JQ, returns the equations on homology
satisfied by a homomorphism between the two corresponding abelian varieties.}

/* Basic invariants */
RR := BaseRing(JP);
gP := #Rows(JP) div 2; gQ := #Rows(JQ) div 2;
n := 4 * gP * gQ;
/* Building a formal matrix corresponding to all possible integral
 * transformations of the lattice */
R := PolynomialRing(RR, n);
vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gQ, 2 * gP, vars);
/* Condition that integral transformation preserve the complex structure */
Comm := Eltseq(M * ChangeRing(JP, R) - ChangeRing(JQ, R) * M);
/* Splitting previous linear equations by formal variable */
return Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] : var in vars ]);

end intrinsic;


intrinsic TangentRepresentation(R::., P::., Q::.) -> .
{Given a homology representation R of a homomorphism and two period matrices P
and Q, returns an analytic representation A of that same homomorphism, so that
A P = Q R.}

CC := BaseRing(P);
P0, s0 := InvertibleSubmatrix(P : IsPeriodMatrix := true);
QR := Q * ChangeRing(R, CC);
QR0 := Submatrix(QR, [ 1..#Rows(QR) ], s0);
/* Invert and return; transposes intervene because of right action */
A := NumericalLeftSolve(P0, QR0);
if Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) gt CC`epscomp then
    error "Error in determining tangent representation";
end if;
return A;

end intrinsic;


intrinsic TangentRepresentation(R::., P::.) -> .
{Given a homology representation R of an endomorphism and a period matrix P,
returns an analytic representation of that same endomorphism, so that A P = Q
R.}

return TangentRepresentation(R, P, P);

end intrinsic;


intrinsic HomologyRepresentation(A::., P::., Q::.) -> .
{Given a complex tangent representation A of a homomorphism and two period
matrices P and Q, returns a homology representation R of that same
homomorphism, so that A P = Q R.}

CC := BaseRing(P);
SplitAP := VerticalSplitMatrix(A * P);
SplitQ := VerticalSplitMatrix(Q);
RRR := NumericalRightSolve(SplitQ, SplitAP);
R := Matrix(Integers(), [ [ Round(cRR) : cRR in Eltseq(row) ] : row in Rows(RRR) ]);
if Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) gt CC`epscomp then
    error "Error in determining homology representation";
end if;
return R;

end intrinsic;


intrinsic HomologyRepresentation(A::., P::.) -> .
{Given a complex tangent representation A of an endomorphism and a period
matrix P, returns a homology representation R of that same endomorphism, so
that A P = P R.}

return HomologyRepresentation(A, P, P);

end intrinsic;


intrinsic GeometricHomomorphismRepresentationPartial(P::., Q::.) -> SeqEnum
{Given period matrices P and Q, this function determines a ZZ-basis of
homomorphisms between the corresponding abelian varieties. These are returned
as pairs of a complex tangent representation A and a homology representation R
for which A P = Q R.}

/* Basic invariants */
gP := #Rows(P); gQ := #Rows(Q);
JP := ComplexStructure(P); JQ := ComplexStructure(Q);

/* Determination of approximate endomorphisms by LLL */
M := RationalHomomorphismEquations(JP, JQ);
Ker := IntegralLeftKernel(M);

/* Deciding which rows to keep */
CC := BaseRing(P); RR := BaseRing(JP); gens := [ ];
for r in Rows(Ker) do
    R := Matrix(Rationals(), 2*gQ, 2*gP, Eltseq(r));
    /* Culling the correct transformations from holomorphy condition */
    Comm := ChangeRing(R, RR) * JP - JQ * ChangeRing(R, RR);
    if &and([ (RR ! Abs(c)) lt RR`epscomp : c in Eltseq(Comm) ]) then
        A := TangentRepresentation(R, P, Q);
        if Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) gt CC`epscomp then
            error "Error in determining homomorphism representations";
        end if;
        Append(~gens, [* A, R *]);
    end if;
end for;
return gens;

end intrinsic;


intrinsic GeometricEndomorphismRepresentationPartial(P::.) -> .
{Given a period matrix P, this function determines a ZZ-basis of endomorphisms
of the corresponding abelian variety. These are returned as pairs of a complex
tangent representation A and a homology representation R for which A P = P R.}

return GeometricHomomorphismRepresentationPartial(P, P);

end intrinsic;


intrinsic GeometricHomomorphismRepresentation(P::., Q::., F::Fld : UpperBound := Infinity()) -> SeqEnum
{Given period matrices P and Q, as well as a field F, this function determines
a ZZ-basis of homomorphisms between the corresponding abelian varieties. These
are returned as triples of an algebraized tangent representation A, a homology
representation R and a complex tangent representation ACC. We have ACC P = Q R,
and via the infinite place of F the matrix A is mapped to ACC.}

gensPart := GeometricHomomorphismRepresentationPartial(P, Q);
gensPol, test := RelativeMinimalPolynomials(gensPart, F : UpperBound := UpperBound);
if not test then
    error "No suitable minimal polynomial found";
end if;
L := RelativeSplittingFieldExtra(gensPol);
gens := [ ];
for gen in gensPart do
    genApp, genHom := Explode(gen);
    genTan, test := AlgebraizeMatrixInRelativeField(genApp, L);
    if not test then
        error "No suitable algebraic element found";
    end if;
    Append(~gens, [* genTan, genHom, genApp *]);
end for;
return gens;

end intrinsic;


intrinsic GeometricEndomorphismRepresentation(P::., F::Fld : UpperBound := Infinity()) -> SeqEnum
{Given a period matrix P and a field F, this function determines a ZZ-basis of
endomorphisms of the corresponding abelian variety. These are returned as
triples of an algebraized tangent representation A, a homology representation R
and a complex tangent representation ACC. We have ACC P = P R, and via the
infinite place of F the matrix A is mapped to ACC.}

return GeometricHomomorphismRepresentation(P, P, F : UpperBound := UpperBound);

end intrinsic;


intrinsic GeometricEndomorphismRepresentationRecognition(gensPart::SeqEnum, L::Fld) -> .
{Given the output of GeometricEndomorphismRepresentationPartial and a field L
with an infinite places, writes the tangent representation of the former as
matrices over the latter. Use when delegating to Sage.}

gens := [ ];
for gen in gensPart do
    genApp, genHom := Explode(gen);
    genTan, test := AlgebraizeMatrixInRelativeField(genApp, L);
    if not test then
        error "No suitable algebraic element found";
    end if;
    Append(~gens, [* genTan, genHom, genApp *]);
end for;
return gens;

end intrinsic;
