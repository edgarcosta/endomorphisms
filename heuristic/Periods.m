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


intrinsic PeriodMatrix(eqsCC::SeqEnum : HaveOldenburg := false) -> AlgMatElt
{Computes a (big) period matrix of the curve defined by the complex polynomials
eqsCC.}

RCC := Parent(eqsCC[1]);
if #GeneratorsSequence(RCC) eq 1 then
    if #eqsCC eq 2 then
        fCC, hCC := Explode(eqsCC);
        gCC := 4*fCC + hCC^2;
    else
        gCC := Explode(eqsCC);
    end if;
    if not HaveOldenburg then
        JCC := AnalyticJacobian(gCC);
        return Transpose(Matrix(BaseRing(gCC), BigPeriodMatrix(JCC)));
    end if;
    return Transpose(Matrix(BaseRing(gCC), PeriodMatrix(gCC)));
elif #GeneratorsSequence(RCC) eq 3 then
    F := Explode(eqsCC);
    SCC<x0,x1,x2> := Parent(F); CC := BaseRing(SCC); RCC<x,y> := PolynomialRing(CC, 2);
    h := hom<SCC -> RCC | [x,y,1]>; f := h(F);
    if not HaveOldenburg then
        error "No functionality for plane curves available";
    end if;
    return Transpose(Matrix(BaseRing(Parent(F)), PeriodMatrix(f)));
else
    error "No functionality for general curves available";
end if;

end intrinsic;
