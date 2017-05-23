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


intrinsic FactorReconstruct(Lat::., K::Fld) -> Crv
{Recovers algebraic expressions for the factors.}

g := #Rows(Transpose(Lat));
if g eq 1 then
    return FactorReconstructG1(Lat, K);
elif g eq 2 then
    return FactorReconstructG2(Lat, K);
else
    error "Reconstruction in genus larger than 2 not implemented yet";
end if;

end intrinsic;


intrinsic FactorReconstructG1(P::., K::Fld) -> .
{Reconstructs elliptic curve factor from a period lattice.}

P := Eltseq(P); CC := Parent(P[1]); RR := RealField(CC);
if Im(P[2]/P[1]) lt 0 then
    P := [ P[2], P[1] ];
end if;
g4CC := 120 * (1/P[1])^4 * ZetaFunction(RR, 4) * Eisenstein(4, P);
g6CC := 280 * (1/P[1])^6 * ZetaFunction(RR, 6) * Eisenstein(6, P);
g4 := AlgebraizeElementInRelativeField(g4CC, K);
g6 := AlgebraizeElementInRelativeField(g6CC, K);
// Division by 16 because of our conventions on period matrices
R<x> := PolynomialRing(K); f := (4*x^3 - g4*x - g6)/16; h := 0;
return HyperellipticCurve(f, h);

end intrinsic;


intrinsic FactorReconstructG2(P::., K::Fld) -> .
{Reconstructs genus 2 factor.}

// Recover small period matrix:
CC := BaseRing(P); RCC<xCC> := PolynomialRing(CC);
g := #Rows(Transpose(P));
Omega1 := Submatrix(P, 1, 1, g, g); Omega2 := Submatrix(P, g + 1, 1, g, g);
PSmall := Omega2 * Omega1^(-1);

rosensCC := RosenhainInvariants(PSmall);
fCC := xCC * (xCC - 1) * &*[ xCC - rosenCC : rosenCC in rosensCC ];
g2sCC := G2Invariants(HyperellipticCurve(fCC));
g2s := [ AlgebraizeElementInRelativeField(g2CC, K) : g2CC in g2sCC ];
X := HyperellipticCurveFromG2Invariants(g2s);
if Type(BaseRing(X)) eq FldRat then
    X := ReducedMinimalWeierstrassModel(X);
end if;
X := TwistDifferentialBasis(X, P);
return X;

end intrinsic;


intrinsic TwistDifferentialBasis(X::., P::.) -> .
{Twists the curve and differential basis as needed.}

F := BaseRing(X);
CC := BaseRing(P); RCC := PolynomialRing(CC);
f, h := HyperellipticPolynomials(X);
fCC := EmbedAtInfinitePlace(f, RCC); hCC := EmbedAtInfinitePlace(h, RCC);
Q := PeriodMatrix([fCC, hCC] : HaveOldenburg:= true);
GeoIsogRep := GeometricIsogenyRepresentationPartial(P, Q);
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
Q := PeriodMatrix([ gCC, 0 ] : HaveOldenburg := true);
GeoIsogRep := GeometricIsogenyRepresentationPartial(P, Q);
A := GeoIsogRep[1][1];

d := AlgebraizeElementInRelativeField(A[1,1]^2, F);
return HyperellipticCurve(d * g);

end intrinsic;


intrinsic FactorDescription(X::Crv) -> List
{Describes factor by recursive list of strings and integers.}

if Type(X) eq CrvHyp then
    return FactorDescriptionHyperelliptic(X);
elif Type(X) eq CrvPln then
    return FactorDescriptionPlane(X);
else
    error "No implementation for general curves yet";
end if;

end intrinsic;


intrinsic FactorDescriptionHyperelliptic(X::CrvHyp) -> List
{Describes factor by recursive list of strings and integers.}

entry1 := "hyp";

// TODO: Fractions may occur. Also, relative case will need slight changes.
F := BaseRing(X);
F_seq := Eltseq(MinimalPolynomial(F.1));
F_seq := [ Rationals() ! coeff : coeff in F_seq ];
// NOTE for relative case:
//F_seq := [ [ Rationals() ! c : c in coeff ] : coeff in F_seq ];
entry2 := F_seq;

f, h := HyperellipticPolynomials(X);
f_seq := Eltseq(f); h_seq := Eltseq(h);
f_seq_seq := [ [ Rationals() ! c : c in Eltseq(coeff) ] : coeff in f_seq ];
h_seq_seq := [ [ Rationals() ! c : c in Eltseq(coeff) ] : coeff in h_seq ];
entry3 := [ f_seq_seq, h_seq_seq ];

L := [* entry1, entry2, entry3 *];
return L;

end intrinsic;


intrinsic FactorDescriptionPlane(X::CrvPln) -> List
{Describes factor by recursive list of strings and integers.}

entry1 := "pln";

// TODO: Fractions may occur. Also, relative case will need slight changes.
F := BaseRing(X);
F_seq := Eltseq(MinimalPolynomial(F.1));
F_seq := [ Rationals() ! coeff : coeff in F_seq ];
// NOTE for relative case:
//F_seq := [ [ Rationals() ! c : c in coeff ] : coeff in F_seq ];
entry2 := F_seq;

f := DefiningPolynomials(X)[1];
mons := Monomials(f);
entry3 := [ ];
for mon in mons do
    coeff := MonomialCoefficient(f, mon);
    coeff_seq := [ Rationals() ! c : c in Eltseq(coeff) ];
    exp := Exponents(mon);
    Append(~entry3, [* coeff_seq, exp *]);
end for;

L := [* entry1, entry2, entry3 *];
return L;

end intrinsic;
