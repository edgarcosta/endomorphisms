/***
 *  Test if given elliptic curve is a factor of the Jacobian
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic TestEllipticFactor(X::Crv, E::Crv) -> BoolElt
{Given a curve X and an elliptic curve E, determines whether E is a factor of
the Jacobian of X. Returns a map from X to E if this is the case.}

if Type(X) eq CrvHyp then
    return TestEllipticFactorHyperelliptic(X, E);
elif Type(X) eq CrvPln then
    return TestEllipticFactorPlane(X, E);
else
    error "Not implemented for general curves yet";
end if;

/* Do this well; there should be no need to force a rational base point, which
 * is only good for endomorphisms at any rate. So perhaps we do not need a
 * distinction. It also helps for debugging purposes. Of course the elliptic
 * factor will have a base point. */
/* The sole remaining complication is the field of definition of the map on
 * differentials itself. The flow should be to determine a field for E, take a
 * compositum (one of the fields being normal) and extending the infinite
 * place, then have the algorithms roll and extend by themselves. This too will
 * be good to write down the API. */
/* Consider triple extension issue at same time */

end intrinsic;


intrinsic TestEllipticCMFactor(X::Crv, D::RngIntElt : prec := 1000) -> BoolElt
{Given a curve X and a discriminant D, determines whether the CM isogeny class
of discriminant D is a factor of the Jacobian of X. Returns a map from X to a
corresponding elliptic curve if this is the case.}

CC := ComplexFieldExtra(prec); RR := RealField(CC);
CCLarge := ComplexFieldExtra(prec + 100);
if D mod 4 ne 0 then
    tau := (Sqrt(CCLarge ! D) + 1)/2;
else
    tau := Sqrt(CCLarge ! D)/2;
end if;
Q := [ CC ! 1, CC ! tau ];

/* TODO: Use this with Shimura or Weber */
g4CC := Eisenstein(4, Q);
g6CC := Eisenstein(6, Q);

QQ := RationalsExtra();
print RelativeMinimalPolynomial(g4CC, QQ);
print RelativeMinimalPolynomial(g6CC, QQ);

return 0;

end intrinsic;
