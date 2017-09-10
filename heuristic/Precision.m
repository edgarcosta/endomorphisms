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
declare attributes FldCom : epscomp, epsLLL, epsinv;
declare attributes FldRe  : epscomp, epsLLL, epsinv;


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Creates a complex field with some extra needed parameters.}

CC := ComplexField(prec);
RR := RealField(CC);
CC`epscomp := RR ! (10^(-prec + 30)); CC`epsLLL  := RR ! (5^(-prec + 2)); CC`epsinv  := RR ! (2^(-prec + 10));
RR`epscomp := RR ! (10^(-prec + 30)); RR`epsLLL  := RR ! (5^(-prec + 2)); RR`epsinv  := RR ! (2^(-prec + 10));
return CC;

end intrinsic;


intrinsic SetEpsComp(CC::FldCom, epscomp::.)
{Modifies epscomp.}

RR := RealField(CC);
CC`epscomp := CC ! epscomp;
RR`epscomp := RR ! epscomp;

end intrinsic;


intrinsic SetEpsLLL(CC::FldCom, epsLLL::.)
{Modifies epsLLL.}

RR := RealField(CC);
CC`epsLLL := CC ! epsLLL;
RR`epsLLL := RR ! epsLLL;

end intrinsic;


intrinsic SetEpsInv(CC::FldCom, epsinv::.)
{Modifies epsinv.}

RR := RealField(CC);
CC`epsinv := CC ! epsinv;
RR`epsinv := RR ! epsinv;

end intrinsic;
