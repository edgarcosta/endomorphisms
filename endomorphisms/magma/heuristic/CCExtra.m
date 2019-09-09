/***
 *  Richer structure for complex number fields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


/* We add parameters for comparison, LLL, and seeing whether a square matrix is
 * invertible. */
declare attributes FldCom : epscomp, epsinv, height_bound;
declare attributes FldRe  : epscomp, epsinv, height_bound;

declare verbose EndoFind, 3;

/* TODO: Some attributes have severe side effects, mainly because there is only one rational field in Magma */


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Returns a complex field of the given precision with the extra attributes
epscomp, epsLLL, epsinv and height_bound.}

CC := ComplexField(prec);
RR := RealField(CC);
CC`epscomp := RR ! (10^(Round(-9.2*prec/10))); CC`epsinv  := RR ! (2^(-prec)); CC`height_bound := RR ! (3^(40 + (prec div 10)));
RR`epscomp := CC`epscomp; RR`epsinv := CC`epsinv; RR`height_bound := CC`height_bound;
return CC;

end intrinsic;


intrinsic ComplexFieldExtra() -> FldCom
{Default ComplexFieldExtra with precision 100.}

return ComplexFieldExtra(100);

end intrinsic;


intrinsic SetEpsComp(CC::FldCom, epscomp::.)
{Modifies the attributes epscomp of CC. This decides whether something is 0.}

RR := RealField(CC);
CC`epscomp := RR ! epscomp;
RR`epscomp := RR ! epscomp;

end intrinsic;


intrinsic SetEpsInv(CC::FldCom, epsinv::.)
{Modifies the attributes epsinv of CC. This decides if a square matrix is invertible via the determinant.}

RR := RealField(CC);
CC`epsinv := RR ! epsinv;
RR`epsinv := RR ! epsinv;

end intrinsic;


intrinsic SetHeightBound(CC::FldCom, height_bound::.)
{Modifies the attributes height_bound of CC. This is the largest coefficient allowed in numerical kernels.}

RR := RealField(CC);
CC`height_bound := RR ! height_bound;
RR`height_bound := RR ! height_bound;

end intrinsic;
