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
declare attributes FldCom : epscomp, epsLLL, epsinv, height_bound;
declare attributes FldRe  : epscomp, epsLLL, epsinv, height_bound;

declare verbose EndoFind, 1;


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Returns a complex field of the given precision with the extra attributes
epscomp, epsLLL, epsinv and height_bound.}

CC := ComplexField(prec);
RR := RealField(CC);
CC`epscomp := RR ! (10^(-prec + 30)); CC`epsLLL  := RR ! (5^(-prec)); CC`epsinv  := RR ! (2^(-prec)); CC`height_bound := RR ! (3^(prec div 2));
//CC`epscomp := RR ! (10^(-prec + 100)); CC`epsLLL  := RR ! (10^(-prec)); CC`epsinv  := RR ! (2^(-prec)); CC`height_bound := RR ! (3^(prec div 2));
//CC`epscomp := RR ! (10^(-prec + 100)); CC`epsLLL  := RR ! (10^(-prec div 2)); CC`epsinv  := RR ! (2^(-prec)); CC`height_bound := RR ! (3^(prec div 2));
RR`epscomp := CC`epscomp; RR`epsLLL := CC`epsLLL; RR`epsinv := CC`epsinv; RR`height_bound := CC`height_bound;
return CC;

end intrinsic;


intrinsic SetEpsComp(CC::FldCom, epscomp::.)
{Modifies the attributes epscomp of CC.}

RR := RealField(CC);
CC`epscomp := RR ! epscomp;
RR`epscomp := RR ! epscomp;

end intrinsic;


intrinsic SetEpsLLL(CC::FldCom, epsLLL::.)
{Modifies the attributes epsLLL of CC.}

RR := RealField(CC);
CC`epsLLL := RR ! epsLLL;
RR`epsLLL := RR ! epsLLL;

end intrinsic;


intrinsic SetEpsInv(CC::FldCom, epsinv::.)
{Modifies the attributes epsinv of CC.}

RR := RealField(CC);
CC`epsinv := RR ! epsinv;
RR`epsinv := RR ! epsinv;

end intrinsic;


intrinsic SetHeightBound(CC::FldCom, height_bound::.)
{Modifies the attributes height_bound of CC.}

RR := RealField(CC);
CC`height_bound := RR ! height_bound;
RR`height_bound := RR ! height_bound;

end intrinsic;
