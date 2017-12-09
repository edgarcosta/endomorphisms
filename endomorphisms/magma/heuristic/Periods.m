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
        gCC := (4*fCC + hCC^2) / 4;
    else
        gCC := Explode(eqsCC);
    end if;
    if not HaveOldenburg then
        JCC := AnalyticJacobian(gCC);
        /* We divide by 2 because we integrate wrt x^i dx / 2y */
        return Transpose(ChangeRing(BigPeriodMatrix(JCC), CC)) / 2;
    end if;
    return Transpose(ChangeRing(PeriodMatrix(gCC : Prec := Precision(CC)), CC)) / 2;
elif #GeneratorsSequence(RCC) eq 3 then
    if not HaveOldenburg then
        error "No functionality for plane curves available";
    end if;
    test, fCC, e := IsSuperelliptic(eqsCC);
    if test then
        P := Transpose(ChangeRing(PeriodMatrix(fCC, e : Prec := Precision(CC)), CC));
        P := SuperellipticCompatibility(P, e);
        return P;
    else
        F := Explode(eqsK);
        X := PlaneCurve(F); f := DefiningEquation(AffinePatch(X, 1));
        return Transpose(ChangeRing(PeriodMatrix(f : Prec := Precision(CC)), CC));
    end if;
else
    error "No functionality for general curves available";
end if;
end intrinsic;


intrinsic IsSuperelliptic(eqs::SeqEnum) -> BoolElt, ., .
{Checks if a curve is superelliptic (in a special form).}

// TODO: Beyond genus 3

R<x,y,z> := Parent(eqs[1]);
if #GeneratorsSequence(R) eq 1 then
    return false, 0, 1;
end if;

F := Explode(eqs);
mons := Monomials(F);
monsy := [ mon : mon in mons | Exponents(mon)[2] ne 0 ];
monsxz := [ mon : mon in mons | Exponents(mon)[2] eq 0 ];
if #monsy ne 1 then
    return false, 0, 1;
end if;

e := Exponents(monsy[1])[2]; C := MonomialCoefficient(F, monsy[1]);
S<t> := PolynomialRing(BaseRing(R));
f := &+[ MonomialCoefficient(F, mon) * t^(Exponents(mon)[1]) : mon in monsxz ];
f := -f/C;
return true, f, e;

end intrinsic;


intrinsic SuperellipticCompatibility(P::., e::RngIntElt) -> .
{Transforms the differentials on a superelliptic curve to compensate for conventions.}

// TODO: Beyond genus 3

Q := Transpose(P);
rowsQ := Rows(Q);
Qtransf := Matrix([ rowsQ[3], rowsQ[1], rowsQ[2] ]);
return Transpose(Qtransf);

end intrinsic;
