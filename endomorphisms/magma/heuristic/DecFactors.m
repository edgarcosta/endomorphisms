/***
 *  Factors from lattices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic FactorReconstruct(P::ModMatFldElt, Q::ModMatFldElt, A::., R::., K::Fld) -> Crv
{Given period matrices P and Q, an analytic and homological representation A
and R of a map between the corresponding Jacobians, and a field K, returns an
algebraic expression for the factor corresponding to Q.}

g := #Rows(Q);
if g eq 1 then
    return FactorReconstructG1(P, Q, A, R, K);
elif g eq 2 then
    return FactorReconstructG2(P, Q, A, R, K);
else
    return 0;
end if;

end intrinsic;


intrinsic FactorReconstructG1(P::ModMatFldElt, Q::ModMatFldElt, A::., R::., K::Fld) -> .
{Given period matrices P and Q, an analytic and homological representation A
and R of a map between the corresponding Jacobians, and a field K, returns an
algebraic expression for the factor corresponding to Q. Assumes that the genus
of the curve corresponding to Q equals 1.}

Q := Eltseq(Q); CC := Parent(Q[1]); RR := RealField(CC);
if Im(Q[2]/Q[1]) lt 0 then
    Q := [ Q[2], Q[1] ];
end if;
g4CC := 120 * (1/Q[1])^4 * ZetaFunction(RR, 4) * Eisenstein(4, Q);
g6CC := 280 * (1/Q[1])^6 * ZetaFunction(RR, 6) * Eisenstein(6, Q);
g4 := AlgebraizeElementInRelativeField(g4CC, K);
g6 := AlgebraizeElementInRelativeField(g6CC, K);
R<x> := PolynomialRing(K); f := (4*x^3 - g4*x - g6)/4; h := 0;
return HyperellipticCurve(f, h);

end intrinsic;


intrinsic FactorReconstructG2(P::ModMatFldElt, Q::ModMatFldElt, A::., R::., K::Fld) -> .
{Given period matrices P and Q, an analytic and homological representation A
and R of a map between the corresponding Jacobians, and a field K, returns an
algebraic expression for the factor corresponding to Q. Assumes that the genus
of the curve corresponding to Q equals 2.}

/* TODO: A lot of improvement needed! */
return 0;
CC := BaseRing(P); RCC<xCC> := PolynomialRing(CC);
gY := #Rows(Transpose(Q));
//return Max([ Abs(c) : c in Eltseq(P*Transpose(A) - ChangeRing(Transpose(R), CC)*Q) ]);
EQ := FindPolarizationBasis(Q)[1];
T := PrincipalBasisRandom(EQ);
Qnew := ChangeRing(T, CC) * Q;
Omega1 := Submatrix(Qnew, 1, 1, gY, gY); Omega2 := Submatrix(Qnew, gY + 1, 1, gY, gY);
tau := Omega2 * Omega1^(-1);
tauIm := Matrix(RealField(CC), [ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
if not IsPositiveDefinite(ChangeRing(tauIm, RealField(Precision(CC) - 10))) then
    tau *:= -1;
end if;
rosensCC := RosenhainInvariants(tau);
fCC := xCC * (xCC - 1) * &*[ xCC - rosenCC : rosenCC in rosensCC ];
g2sCC := G2Invariants(HyperellipticCurve(fCC));
return g2sCC;

g2s := [ AlgebraizeElementInRelativeField(g2CC, K) : g2CC in g2sCC ];
Y := HyperellipticCurveFromG2Invariants(g2s);
if Type(BaseRing(Y)) eq FldRat then
    Y := ReducedMinimalWeierstrassModel(Y);
end if;
Y := TwistDifferentialBasis(Y, Qnew);
return Y;

end intrinsic;


intrinsic TwistDifferentialBasis(X::., P::.) -> .
{Given a period matrix P and a curve X that represents P over CC, twists X so
that it also represents P over the base field.}
/* TODO: not for now, and only applies on generic stratum */

F := BaseRing(X);
CC := BaseRing(P); RCC := PolynomialRing(CC);
f, h := HyperellipticPolynomials(X);
fCC := EmbedAtInfinitePlace(f, RCC); hCC := EmbedAtInfinitePlace(h, RCC);
Q := PeriodMatrix([fCC, hCC] : MolinNeurohr:= true);
GeoIsogRep := GeometricHomomorphismRepresentationCC(P, Q);
A := GeoIsogRep[1][1];

pols := [ RelativeMinimalPolynomial(c, F) : c in Eltseq(A) ];
sqrtCC := CC ! 1; d := 1;
for pol in pols do
    if Degree(pol) eq 2 then
        polCC := EmbedAtInfinitePlace(pol, RCC);
        a := Coefficient(pol, 2); b := Coefficient(pol, 1); c := Coefficient(pol, 0);
        aCC := Coefficient(polCC, 2); bCC := Coefficient(polCC, 1); cCC := Coefficient(polCC, 0);
        sqrtCC := Sqrt(bCC^2 - 4*aCC*cCC); d := b^2 - 4*a*c;
    end if;
end for;

WA := AlgebraizeMatrixInRelativeField(A / sqrtCC, F);
g := 4*f + h^2; R<x> := Parent(g);
num := WA[2,2]*x + WA[1,2]; den := WA[2,1]*x + WA[1,1];
g := R ! (den^6 * Evaluate(g, num/den));
gCC := EmbedAtInfinitePlace(g, RCC);
Q := PeriodMatrix([ gCC, 0 ] : MolinNeurohr := true);
GeoIsogRep := GeometricHomomorphismRepresentationCC(P, Q);
A := GeoIsogRep[1][1];

d := AlgebraizeElementInRelativeField(A[1,1]^2, F);
return HyperellipticCurve(d * g);

end intrinsic;


intrinsic FactorDescription(X::Crv, F::Fld) -> List
{Returns a string description of the curve X over the field F.}

if Type(X) eq CrvHyp then
    return FactorDescriptionHyperelliptic(X, F);
elif Type(X) eq CrvPln then
    return FactorDescriptionPlane(X, F);
else
    error "No implementation for general curves yet";
end if;

end intrinsic;


intrinsic FactorDescriptionHyperelliptic(X::CrvHyp, F::Fld) -> List
{Returns a string description of the curve X over the field F.}

if Genus(X) eq 1 then
    desc := "ell";

    K := BaseRing(X); F := BaseRing(K);
    K_seq := FieldDescription(K, F);
    field := K_seq;

    f, h := HyperellipticPolynomials(X);
    f_seq := Eltseq(f); h_seq := Eltseq(h);
    f_seq_seq := [ ElementDescription(coeff, F) : coeff in f_seq ];
    h_seq_seq := [ ElementDescription(coeff, F) : coeff in h_seq ];
    coeffs := [ f_seq_seq, h_seq_seq ];

else
    desc := "hyp";
    field := [ ];
    coeffs := [ ];
end if;

return [* desc, field, coeffs *];

end intrinsic;


intrinsic FactorDescriptionPlane(X::CrvPln, F::Fld) -> List
{Returns a string description of the curve X over the field F.}

desc := "pln";

K := BaseRing(X); F := BaseRing(K);

K_seq := FieldDescription(K, F);
field := K_seq;

f := DefiningPolynomials(X)[1];
mons := Monomials(f);
coeffsexps := [ ];
for mon in mons do
    coeff := MonomialCoefficient(f, mon);
    coeff_seq := ElementDescription(coeff, F);
    exp := Exponents(mon);
    Append(~coeffsexps, [* coeff_seq, exp *]);
end for;

return [* desc, field, coeffsexps *];

end intrinsic;
