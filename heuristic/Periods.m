/***
 *  Determining period matrices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


/* Enable Oldenburg if you have access to the relevant code by Pascal Molin,
 * Christian Neurohr et al. */


intrinsic PeriodMatrix(eqsCC::SeqEnum, eqsK::SeqEnum : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the curve defined by the complex polynomials
eqsCC.}

RCC := Parent(eqsCC[1]); CC := BaseRing(RCC);
if #GeneratorsSequence(RCC) eq 1 then
    if #eqsCC eq 2 then
        fCC, hCC := Explode(eqsCC);
        gCC := 4*fCC + hCC^2;
    else
        gCC := Explode(eqsCC);
    end if;
    if not HaveOldenburg then
        JCC := AnalyticJacobian(gCC);
        return Transpose(Matrix(CC, BigPeriodMatrix(JCC)));
    end if;
    return Transpose(Matrix(CC, PeriodMatrix(gCC : Prec := Precision(CC))));
elif #GeneratorsSequence(RCC) eq 3 then
    if not HaveOldenburg then
        error "No functionality for plane curves available";
    end if;
    F := Explode(eqsK);
    S<x0,x1,x2> := Parent(F); K := BaseRing(S); R<x,y> := PolynomialRing(K, 2);
    h := hom<S -> R | [x,y,1]>; f := h(F);
    return Transpose(Matrix(CC, PeriodMatrix(f : Prec := Precision(CC))));
    /*
    FCC := Explode(eqsCC);
    SCC<x0,x1,x2> := Parent(FCC); CC := BaseRing(SCC); RCC<x,y> := PolynomialRing(CC, 2);
    h := hom<SCC -> RCC | [x,y,1]>; fCC := h(FCC);
    return Transpose(Matrix(CC, PeriodMatrix(fCC)));
    */
else
    error "No functionality for general curves available";
end if;

end intrinsic;
