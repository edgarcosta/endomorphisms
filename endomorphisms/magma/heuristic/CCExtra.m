/***
 *  Richer structure for complex number fields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


/* We add parameters for comparison, LLL, and seeing whether a square matrix is
 * invertible. */
declare attributes FldCom : epscomp, epsinv, height_bound, prec_algdep;
declare attributes FldRe  : epscomp, epsinv, height_bound, prec_algdep;

declare verbose EndoFind, 3;

// TODO: We may want to have three different fields: default field, period
// matrix, iota, in that order of precision.
// Right now this is only asking for work though.


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Returns a complex field of the given precision with the extra attributes
epscomp, epsLLL, epsinv and height_bound.}

CC := ComplexField(prec);
RR := RealField(CC);
if prec lt 200 then
    CC`epscomp := RR ! (10^(-prec + 10)); CC`epsinv  := RR ! 10^(-6); CC`height_bound := 10^6; CC`prec_algdep := prec - 5;
    RR`epscomp := CC`epscomp; RR`epsinv := CC`epsinv; RR`height_bound := CC`height_bound; RR`prec_algdep := CC`prec_algdep;
    return CC;
end if;
CC`epscomp := RR ! (10^(-Round(9.5*prec/10))); CC`epsinv  := RR ! (2^(-prec)); CC`height_bound := RR ! (3^(30 + (prec div 10))); CC`prec_algdep := Round(8*prec/10);
RR`epscomp := CC`epscomp; RR`epsinv := CC`epsinv; RR`height_bound := CC`height_bound; RR`prec_algdep := CC`prec_algdep;
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


intrinsic AlmostEqual(approx::., alg::. : iota:=0) -> BoolElt
{ return if approx is floating point approximation of alg. The second return value is the relative error }
    R := Parent(approx);
    require assigned R`epscomp: "the parent of the first argument must be a {Real,Complex}FieldExtra";
    Ralg := Parent(alg);
    if Type(Ralg) in [FldCom, FldRe] then // the alg is also a floating point approximation
        eps := Minimum(R`epscomp, Ralg`epscomp);
        maxnorm := Maximum(Abs(approx), Abs(alg));
        normalizer := maxnorm lt eps select 1 else maxnorm;
        algapprox := alg;
    else
        eps := Minimum(R`epscomp, Ralg`CC`epscomp);
        algapprox := EmbedExtra(alg : iota:=iota);
        normalizer := alg eq 0 select 1 else Abs(algapprox);
    end if;
    err := Abs(approx - algapprox)/normalizer;
    return err lt eps, err;
end intrinsic;

intrinsic AlmostEqualMatrix(approx::., alg::. : iota:=0) -> BoolElt
{ return if approx is an approximation of r }
    require assigned BaseRing(approx)`epscomp: "the base ring of th first argument must be RealFieldExtra or ComplexFieldExtra";
    require forall{ f : f in [Nrows, Ncols] | f(approx) eq f(alg) } : "the matrices must have the same dimension";
    approxseq := Eltseq(approx);
    algseq := Eltseq(alg);
    err := 0;
    // we don't use forall to track the error
    res := [<b,e> where b, e := AlmostEqual(approxseq[i], algseq[i] : iota:=iota) : i in [1..#algseq]];
    return &and[elt[1] : elt in res], Maximum([elt[2] : elt in res]);
end intrinsic;

//SetDebugOnError(true);R<x>:=PolynomialRing(Rationals()); C:=HyperellipticCurve(R![96,0,754,0,-118,0,-48]); HeuristicDecomposition(C);
