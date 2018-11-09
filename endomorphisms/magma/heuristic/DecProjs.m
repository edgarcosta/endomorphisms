/***
 *  Decomposing Jacobians
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


forward MatricesFromIdempotent;
forward NonCentralIdempotentsStepOne;
forward NonCentralIdempotentsStepTwo;
forward FindIdempotentsStupid;
forward IsTrueIdempotent;
/* TODO: Too much information is recalculated: this should be reproved by true abelian functionality if possible */
/* TODO: Endomorphism calculation of factors should proceed by using maps to factor and its dual. */


intrinsic IsotypicalIdempotents(P::., EndoRep::.) -> .
{Returns factors of the Jacobian and appropriate spanning set of idempotents over CC.}

EndoAlg, EndoDesc := EndomorphismStructure(EndoRep);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
C := EndoAlg[1]; idemsC := CentralIdempotents(C);
return [ MatricesFromIdempotent(idemC, EndoData) : idemC in idemsC ];

end intrinsic;


function MatricesFromIdempotent(idem, EndoData)
// Returns the matrix representations corresponding to the idempotents in
// idems, using the endomorphism structure EndoData.

EndoRep, EndoAlg, EndoDesc := Explode(EndoData);
As := [ gen[1] : gen in EndoRep ]; Rs := [ gen[2] : gen in EndoRep ];
C, GensC := Explode(EndoAlg);

idemC := [ Rationals() ! c : c in Eltseq(idem) ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
test, idem := MatrixInBasis(idemC, GensC);
assert test;

idem := Eltseq(idem);
idemA := &+[ idem[i] * As[i] : i in [1..#As] ];
idemR := &+[ idem[i] * Rs[i] : i in [1..#Rs] ];
return [* idemA, idemR *];

end function;


intrinsic ComponentWithProjection(P::., idem::.) -> .
{Returns isotypical component B^d of the Jacobian and projection to it.}

K := BaseRing(idem[1]);
R := idem[2];
R := R*Denominator(R);
ACC := TangentRepresentation(R, P, P);
Q, h := ImgProj([* ACC, R *], P, P);
BCC := h[1]; S := h[2];
/* Recalculation */
test, B := AlgebraizeMatrix(BCC, K);
assert test;
return Q, [* B, S *];

end intrinsic;


intrinsic ComponentWithInclusion(P::., idem::.) -> .
{Returns isotypical component B^d of the Jacobian and inclusion from it.}

K := BaseRing(idem[1]);
R := 1 - idem[2];
R := R*Denominator(R);
ACC := TangentRepresentation(R, P, P);
Q, h := Ker0([* ACC, R *], P, P);
BCC := h[1]; S := h[2];
/* Recalculation */
test, B := AlgebraizeMatrix(BCC, K);
assert test;
return Q, [* B, S *];

end intrinsic;


intrinsic IsotypicalComponentsWithProjections(P::., EndoRep::.) -> .
{Returns isotypical components B^d of the Jacobian and projections to them.}

idems := IsotypicalIdempotents(P, EndoRep);
comps := [* *];
for idem in idems do
    Q, h := ComponentWithProjection(P, idem); 
    Append(~comps, [* Q, h *]);
end for;
return comps;

end intrinsic;


intrinsic IsotypicalComponentsWithInclusions(P::., EndoRep::.) -> .
{Returns isotypical components B^d of the Jacobian and inclusions to it.}

idems := IsotypicalIdempotents(P, EndoRep);
comps := [* *];
for idem in idems do
    Q, h := ComponentWithInclusion(P, idem); 
    Append(~comps, [* Q, h *]);
end for;
return comps;

end intrinsic;


intrinsic IsotypicalComponentFoD(Q::., h::.) -> .
{Returns field of definition.}

K := BaseRing(h[1]);
return SubfieldExtra(K, Eltseq(h[1]));

end intrinsic;


intrinsic IsotypicalComponentSplittingIdempotents(Q::., h::.) -> .
{Returns further idempotents. Parent is smallest field over which splitting occurs.}
/* Note that in this case the field gets determined along with the idempotent */

L := BaseRing(h[1]); K := IsotypicalComponentFoD(Q, h);
GeoEndoRepCC := GeometricEndomorphismRepresentationCC(Q);
GeoEndoRep := [ ];
for tupCC in GeoEndoRepCC do
    test, A := AlgebraizeMatrix(tupCC[1], K);
    R := tupCC[2];
    Append(~GeoEndoRep, [* A, R *]);
end for;
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
idems_geo := NonCentralIdempotents(GeoEndoData);

/* Find automorphism group (over QQ for now) */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
AutL := [* Gp, Gf, Gphi *];
H := FixedGroupExtra(L, K, AutL);

S := Subgroups(H); Js := [ rec`subgroup : rec in S ];
Js := [ J : J in Js | #J ne 1 ];
if #Js eq 0 then
    return idems_geo;
end if;

/* Run through lattice */
for J in Reverse(Js) do
    GalM := [* Generators(J), Gphi *];
    EndoData := EndomorphismData(GeoEndoRep, GalM);
    idems := NonCentralIdempotents(EndoData);
    if #idems eq #idems_geo then
        return idems;
    end if;
end for;
    
end intrinsic;


intrinsic RootsOfIsotypicalComponentWithInclusions(Q::., h::.) -> .
{Returns component along with maps. Parent is smallest field over which splitting occurs.}

idems := IsotypicalComponentSplittingIdempotents(Q, h);
L := BaseRing(h[1]); K := BaseRing(idems[1][1]); hKL := CanonicalInclusionMap(K, L);
A, R := Explode(h); A := CoerceToSubfieldMatrix(A, L, K, hKL);
comps := [* *];
for idem in idems do
    Qroot, hnew := ComponentWithInclusion(Q, idem); 
    B, S := Explode(hnew);
    hcomp := [* A*B, R*S *];
    Append(~comps, [* Qroot, hcomp *]);
end for;
return comps;

end intrinsic;


intrinsic NonCentralIdempotents(EndoData::.) -> .
{Returns further decomposition.}

test, idems := NonCentralIdempotentsStepOne(EndoData);
if test then
    return [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
end if;
test, idems := NonCentralIdempotentsStepTwo(EndoData);
if test then
    return [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
end if;
print "All known cases in NonCentralIdempotents fell through, returning unit element";
return [ MatricesFromIdempotent(EndoData[2][1] ! 1, EndoData) : idem in idems ];

end intrinsic;


function NonCentralIdempotentsStepOne(EndoData)

C := EndoData[2][1];
E1, f1 := AlgebraOverCenter(C);
/* This seems a bit heavy-handed */
//F := ClearDenominator(BaseRing(E1));
//if Type(F) eq FldNum then
//    F := ClearDenominator(F);
//    F := ImproveField(F);
//end if;
//E2, f2 := ChangeRing(E1, F);
E2 := E1;

/* Fields have no idempotents */
if IsCommutative(E2) then
    return true, [ C ! 1 ];
end if;

test_dim, d := IsSquare(Dimension(E2));
/* Right now can only deal with algebras of dimension 4 */
if d eq 2 then
    test_quat, Q, f3 := IsQuaternionAlgebra(E2);
    test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
    //f := f1 * f2 * f3 * f4;
    f := f1 * f3 * f4;
    invf := Inverse(f);
    return true, [ C ! invf(M ! [1,0,0,0]), C ! invf(M ! [0,0,0,1]) ];
end if;
return false, [ ];

end function;


function NonCentralIdempotentsAlgebraStepTwo(EndoData);
/* M3 over QQ or CM, M4 over ZZ or CM, M2 over quaternion algebra */

C := EndoData[2][1]; d := Dimension(C);
if d in { 9, 18 } then
    return true, FindIdempotentsStupid(EndoData, 1, 3);
elif d in { 32 } then
    return true, FindIdempotentsStupid(EndoData, 1, 4);
elif d in { 16 } then
    disc := Abs(Discriminant(MaximalOrder(C)));
    if disc eq 1 then
        return true, FindIdempotentsStupid(EndoData, 1, 4);
    else
        return true, FindIdempotentsStupid(EndoData, 2, 2);
    end if;
end if;
return false, [ ];

end function;


function FindIdempotentsStupid(EndoData, dim, goal)

C := EndoData[2][1]; d := Dimension(C);
Bmin := 0; Bmax := 0;
counter := 0;
idems := [ ];

while true do
    if IsEven(counter) then
        Bmax +:= 1;
    else
        Bmin -:= 1;
    end if;
    counter +:= 1;
    D := [ Bmin..Bmax ];

    CP := CartesianPower(D, d);
    for tup in CP do
        c := C ! tup[1];
        v := Eltseq(c); v2 := Eltseq(c^2);
        _, lambda := IsMultiple(v2, v);
        idem := (1/lambda)*c;
        if IsTrueIdempotent(EndoData, dim, idem) then
            test := MatrixInBasis(idem, idems);
            if not test then
                Append(~idems, idem);
            end if;
        end if;
        if #idems eq goal then
            return idems;
        end if;
    end for;
end while;

end function;


function IsTrueIdempotent(EndoData, dim, idem)

A, R := Explode(MatricesFromIdempotent(idem, EndoData));
if Dimension(1 - R) eq 2*dim then
    return true;
end if;
return false;

end function;


intrinsic RootsOfIsotypicalComponentWithProjections(Q::., h::.) -> .
{Returns component along with maps. Parent is smallest field over which splitting occurs.}

idems := IsotypicalComponentSplittingIdempotents(Q, h);
L := BaseRing(h[1]); K := BaseRing(idems[1][1]); hKL := CanonicalInclusionMap(K, L);
A, R := Explode(h); A := CoerceToSubfieldMatrix(A, L, K, hKL);
comps := [* *];
for idem in idems do
    Qroot, hnew := ComponentWithProjection(Q, idem); 
    B, S := Explode(hnew);
    hcomp := [* B*A, S*R *];
    Append(~comps, [* Qroot, hcomp *]);
end for;
return comps;

end intrinsic;


intrinsic ReconstructCurveFromRoot(root::.) -> .
{Curve reconstruction with extension if needed.}

Qroot, hcomp := Explode(root); K := BaseRing(hcomp[1]);
g := #Rows(Qroot);
if g eq 1 then
    return ReconstructCurveG1(Qroot, K);
elif g eq 2 then
    /* TODO: First take isogeny! */
    return ReconstructCurveG2(Qroot, K);
end if;
error "Reconstruction for genus larger than 2 not yet implemented";

end intrinsic;
