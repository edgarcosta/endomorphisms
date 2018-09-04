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


/* Enable MolinNeurohr to use the new code Pascal Molin and Christian Neurohr,
 * which is highly recommended */

intrinsic PeriodMatrix(eqsCC::SeqEnum, eqsK::SeqEnum : MolinNeurohr := true) -> ModMatFldElt
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
        /* We divide by 2 because we integrate with respect to the canonical
         * differential x^i dx / 2y
         * (MN use x^i dx) */
        return ChangeRing(BigPeriodMatrix(JCC), CC) / 2;
    else
        X := SE_Curve(gCC, 2 : Prec := Precision(CC));
        return ChangeRing(X`BigPeriodMatrix, CC) / 2;
    end if;

elif #GeneratorsSequence(RCC) eq 3 then
    if not MolinNeurohr then
        error "No functionality for plane curves available";
    end if;
    test, fCC, e := IsSuperellipticEquation(eqsCC);
    if false then
        X := SE_Curve(fCC, 3 : Prec := Precision(CC));
        P := X`BigPeriodMatrix;
        return SuperellipticCompatibility(P, e);
    else
        /* Note: only polynomials over QQ for now */
        F := Explode(eqsK);
        X := PlaneCurve(F); f := DefiningEquation(AffinePatch(X, 1));
        try 
            /* TODO: Add this when it comes */
            //return ChangeRing(BigPeriodMatrix(RiemannSurface(f : Prec := Precision(CC))), CC);
            return ChangeRing(RS_BigPeriodMatrix(f : Prec := Precision(CC)), CC);
            return 1/(1 - 1);
        catch err
            error "No functionality for plane curves available";
        end try;
    end if;

else
    error "No functionality for general curves available";
end if;
end intrinsic;


intrinsic PeriodMatrix(X::Crv : prec := 300, MolinNeurohr := true) -> ModMatFldElt
{Returns the period matrix of the curve defined by the complex polynomials
eqsCC.}

eqsCC := EmbedCurveEquations(X, prec);
eqsF := DefiningEquations(X);
return PeriodMatrix(eqsCC, eqsF : MolinNeurohr := MolinNeurohr);

end intrinsic;


intrinsic IsSuperellipticEquation(eqs::SeqEnum) -> BoolElt, ., .
{Returns whether the plane curve defined by eqs is of the form y^e z^* = f (x,
z). If so, return the inhomogenous form of f along with e.}

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

e := Exponents(monsy[1])[2];
S<t> := PolynomialRing(BaseRing(R));
f := &+[ MonomialCoefficient(F, mon) * t^(Exponents(mon)[1]) : mon in monsxz ];
C := MonomialCoefficient(F, monsy[1]);
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
