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


/* Enable MolinNeurohr if you have access to the relevant code by Pascal Molin and
 * Christian Neurohr */

intrinsic PeriodMatrix(eqsCC::SeqEnum, eqsK::SeqEnum : MolinNeurohr := false) -> AlgMatElt
{Returns the period matrix of the curve defined by the complex polynomials
eqsCC.}

RCC := Parent(eqsCC[1]); CC := BaseRing(RCC);
if #GeneratorsSequence(RCC) eq 1 then
    if #eqsCC eq 2 then
        fCC, hCC := Explode(eqsCC);
        gCC := (4*fCC + hCC^2) / 4;
    else
        gCC := Explode(eqsCC);
    end if;
    if not MolinNeurohr then
        JCC := AnalyticJacobian(gCC);
        /* We divide by 2 because we integrate wrt x^i dx / 2y */
        return ChangeRing(BigPeriodMatrix(JCC), CC) / 2;
    end if;
    X := SE_Curve(gCC, 2 : Prec := Precision(CC));
    return ChangeRing(X`BigPeriodMatrix, CC) / 2;
    /* Alternative version: */
    //return ChangeRing(PeriodMatrix(gCC : Prec := Precision(CC)), CC) / 2;
elif #GeneratorsSequence(RCC) eq 3 then
    if not MolinNeurohr then
        error "No functionality for plane curves available";
    end if;
    test, fCC, e := IsSuperelliptic(eqsCC);
    if test then
        P := ChangeRing(PeriodMatrix(fCC, e : Prec := Precision(CC)), CC);
        return SuperellipticCompatibility(P, e);
    else
        F := Explode(eqsK);
        X := PlaneCurve(F); f := DefiningEquation(AffinePatch(X, 1));
        return ChangeRing(PeriodMatrix(f : Prec := Precision(CC)), CC);
    end if;
else
    error "No functionality for general curves available";
end if;
end intrinsic;


intrinsic IsSuperelliptic(eqs::SeqEnum) -> BoolElt, ., .
{Returns whether the plane curve defined by eqs is of the form y^e z^* = f (x,
z). If so, return the inhomogenous form of f along with e.}
// TODO: Deal with this beyond genus 3 by making superelliptic curves a class
// of their own.

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
// TODO: Generalize this to apply beyond genus 3. This is a matter of fixing a
// base of differentials. But actually superelliptic curves should be treated
// as a class of their own.

rowsP := Rows(P);
if #rowsP eq 3 then
    return Matrix([ rowsP[3], rowsP[1], rowsP[2] ]);
end if;
error "Need g = 3 for now";

end intrinsic;
