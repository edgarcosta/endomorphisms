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
/* TODO: Endomorphism calculation of factors should instead work by using maps to factor and its dual. */
/* TODO: Make choice between projection and inclusion a flag */


intrinsic IsotypicalIdempotents(P::., GeoEndoRep::.) -> .
{Returns factors of the Jacobian and appropriate spanning set of idempotents over CC.}

GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
C := GeoEndoAlg[1]; idemsC := CentralIdempotents(C);
return [ MatricesFromIdempotent(idemC, GeoEndoData) : idemC in idemsC ];

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


intrinsic ComponentWithProjection(P::., idem::. : CoerceToBase := true) -> .
{Returns projection from Jacobian P corresponding to idempotent idem.}

L := BaseRing(idem[1]);
R := idem[2]; //R := R*Denominator(R);
ACC := TangentRepresentation(R, P, P);
Q, h := ImgProj([* ACC, R *], P, P);
BCC := h[1]; S := h[2];
/* Recalculation */
test, B := AlgebraizeMatrix(BCC, L);
assert test;
K, hKL := SubfieldExtra(L, Eltseq(B));
if CoerceToBase then
    B := CoerceToSubfieldMatrix(B, L, K, hKL);
end if;
incdata := [* L, K, hKL *];
return Q, [* B, S *], incdata;

end intrinsic;


intrinsic ComponentWithInclusion(P::., idem::. : CoerceToBase := true) -> .
{Returns inclusion into Jacobian P corresponding to idempotent idem.}

L := BaseRing(idem[1]);
R := 1 - idem[2]; //R := R*Denominator(R);
ACC := TangentRepresentation(R, P, P);
Q, h := Ker0([* ACC, R *], P, P);
BCC := h[1]; S := h[2];
/* Recalculation */
test, B := AlgebraizeMatrix(BCC, L);
assert test;
if not CoerceToBase then
    return Q, [* B, S *], 0;
end if;
K, hKL := SubfieldExtra(L, Eltseq(B));
if CoerceToBase then
    B := CoerceToSubfieldMatrix(B, L, K, hKL);
end if;
incdata := [* L, K, hKL *];
return Q, [* B, S *], incdata;

end intrinsic;


intrinsic IsotypicalComponentsWithProjections(P::., EndoRep::. : CoerceToBase := true) -> .
{Returns isotypical components Q = B^d of the Jacobian and projections from P to these Q.}

idems := IsotypicalIdempotents(P, EndoRep);
comps := [ ];
for idem in idems do
    Q, h, incdata := ComponentWithProjection(P, idem : CoerceToBase := CoerceToBase);
    Append(~comps, [* Q, h, incdata *]);
end for;
return comps;

end intrinsic;


intrinsic IsotypicalComponentsWithInclusions(P::., EndoRep::. : CoerceToBase := true) -> .
{Returns isotypical components Q = B^d of the Jacobian and inclusions from these Q to P.}

idems := IsotypicalIdempotents(P, EndoRep);
comps := [ ];
for idem in idems do
    Q, h, incdata := ComponentWithInclusion(P, idem : CoerceToBase := CoerceToBase);
    Append(~comps, [* Q, h, incdata *]);
end for;
return comps;

end intrinsic;


intrinsic SplittingIdempotents(Q::., h::., incdata::.) -> .
{Returns further idempotents over the smallest field where the isotypical component splits as far as possible.}

L, K, hKL := Explode(incdata);
/* Recalculate endomorphism algebra over known field (as mentioned above, this is stupid) */
GeoEndoRepCC := GeometricEndomorphismRepresentationCC(Q);
GeoEndoRep := [ ];
for tupCC in GeoEndoRepCC do
    test, A := AlgebraizeMatrix(tupCC[1], K); R := tupCC[2];
    Append(~GeoEndoRep, [* A, R *]);
end for;
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
idems_geo := NonCentralIdempotents(GeoEndoData);

/* Find automorphism group (over QQ for now) */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
AutL := [* Gp, Gf, Gphi *];
H := FixedGroupExtra(L, K, hKL, AutL);

/* Find subgroups corresponding to fields between L and K */
S := Subgroups(H); Js := [ rec`subgroup : rec in S ];
Js := [ J : J in Js | #J ne 1 ];
if #Js eq 0 then
    return idems_geo, [* L, L, CanonicalInclusionMap(L, L) *];
end if;

/* Run through lattice and use smallest field where decomposition occurs */
for J in Reverse(Js) do
    gensM := Generators(J); GalM := [* gensM, Gphi *];
    EndoData := EndomorphismData(GeoEndoRep, GalM);
    idems := NonCentralIdempotents(EndoData);
    if #idems eq #idems_geo then
        M, hML := FixedFieldExtra(L, [ Gphi(genM) : genM in gensM ]);
        return idems, [* L, M, hML *];
    end if;
end for;

end intrinsic;


intrinsic NonCentralIdempotents(EndoData::.) -> .
{Returns further decomposition.}

/* Assert small dimension */
g := #Rows(EndoData[1][1][1]);
assert g le 4;

/* Try successively more complicated methods */
test1, idems := NonCentralIdempotentsStepOne(EndoData);
if test1 then
    return [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
end if;
test2, idems := NonCentralIdempotentsStepTwo(EndoData);
if test2 then
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

/* Central algebras of dimension 4 */
test_dim, d := IsSquare(Dimension(E2));
if d eq 2 then
    test_quat, Q, f3 := IsQuaternionAlgebra(E2);
    test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
    if not test_mat then
        return true, [ C ! 1 ];
    end if;
    //f := f1 * f2 * f3 * f4;
    f := f1 * f3 * f4;
    invf := Inverse(f);
    return true, [ C ! invf(M ! [1,0,0,0]), C ! invf(M ! [0,0,0,1]) ];
end if;
return false, [ ];

end function;


function NonCentralIdempotentsAlgebraStepTwo(EndoData);
/* Catches M3 over QQ or CM, M4 over ZZ or CM, M2 over quaternion algebra */

C := EndoData[2][1]; d := Dimension(C);
dim := 0; goal := 0;
if d in { 9, 18 } then
    dim := 1; goal := 3;
elif d in { 32 } then
    dim := 1; goal := 4;
elif d in { 16 } then
    disc := Abs(Discriminant(MaximalOrder(C)));
    if disc eq 1 then
        dim := 1; goal := 4;
    else
        dim := 2; goal := 2;
    end if;
end if;
if dim ne 0 then
    return true, FindIdempotentsStupid(EndoData, dim, goal);
end if;
return false, [ ];

end function;


function FindIdempotentsStupid(EndoData, dim, goal)

C := EndoData[2][1]; d := Dimension(C);
Bmin := 0; Bmax := 0; counter := 0; idems := [ ];

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
        test_mult, lambda := IsMultiple(v2, v);
        if test_mult then
            idem := (1/lambda)*c;
            if IsTrueIdempotent(EndoData, dim, idem) then
                test_seen := MatrixInBasis(idem, idems);
                if not test_seen then
                    Append(~idems, idem);
                end if;
            end if;
            if #idems eq goal then
                return idems;
            end if;
        end if;
    end for;
end while;

end function;


function IsTrueIdempotent(EndoData, dim, idem)
/* Only accept idempotents whose homology representation has rank dim */

A, R := Explode(MatricesFromIdempotent(idem, EndoData));
if Dimension(1 - R) eq 2*dim then
    return true;
end if;
return false;

end function;


intrinsic RootsOfIsotypicalComponentWithProjections(Q::., h::., incdata::.) -> .
{Returns components along with projection maps over smallest possible field.}
/* We emphatically do not want a single component and multiple maps to it here,
 * because that causes us to miss factors that become isogenous only later */

idems, incdataroot := SplittingIdempotents(Q, h, incdata);
comps := [ ];
for idem in idems do
    /* No need to coerce since we are already over smallest possible field
     * Still, it might be interesting to know what the morphism generates */
    Qroot, hroot, _ := ComponentWithProjection(Q, idem : CoerceToBase := false);
    Append(~comps, [* Qroot, hroot, incdataroot *]);
end for;
return comps;

end intrinsic;


intrinsic RootsOfIsotypicalComponentWithInclusions(Q::., h::., incdata::.) -> .
{Returns components along with inclusion maps over smallest possible field.}
/* We emphatically do not want a single component and multiple maps to it here,
 * because that causes us to miss factors that become isogenous only later */

idems, incdataroot := SplittingIdempotents(Q, h, incdata);
comps := [ ];
for idem in idems do
    /* No need to coerce since we are already over smallest possible field
     * Still, it might be interesting to know what the morphism generates */
    Qroot, hroot, _ := ComponentWithInclusion(Q, idem : CoerceToBase := false);
    Append(~comps, [* Qroot, hroot, incdataroot *]);
end for;
return comps;

end intrinsic;


intrinsic SplitComponentsWithProjections(P::., GeoEndoRep::.) -> .
{Returns maximal possible splitting of the Jacobian P over the smallest field over which this occurs, plus corresponding projections.}

L := BaseRing(GeoEndoRep[1][1]);
comps := [ ];
comps_iso := IsotypicalComponentsWithProjections(P, GeoEndoRep : CoerceToBase := false);
for comp_iso in comps_iso do
    Q, h, incdata := Explode(comp_iso);
    A, R := Explode(h);
    for comp_root in RootsOfIsotypicalComponentWithProjections(Q, h, incdata) do
        Qroot, hroot, incdataroot := Explode(comp_root);
        print incdataroot;
        L, K, hKL := Explode(incdataroot);
        A0 := CoerceToSubfieldMatrix(A, L, K, hKL);
        Aroot, Rroot := Explode(hroot);
        hcomp := [* Aroot*A0, Rroot*R *];
        Append(~comps, [* Qroot, hcomp, incdataroot *]);
    end for;
end for;
return comps;

end intrinsic;


intrinsic SplitComponentsWithInclusions(P::., GeoEndoRep::.) -> .
{Returns maximal possible splitting of the Jacobian P over the smallest field over which this occurs, plus corresponding inclusions.}

L := BaseRing(GeoEndoRep[1][1]);
comps := [ ];
comps_iso := IsotypicalComponentsWithInclusions(P, GeoEndoRep : CoerceToBase := false);
for comp_iso in comps_iso do
    Q, h, incdata := Explode(comp_iso);
    A, R := Explode(h);
    for comp_root in RootsOfIsotypicalComponentWithProjections(Q, h, incdata) do
        Qroot, hroot, incdataroot := Explode(comp_root);
        L, K, hKL := Explode(incdataroot);
        A0 := CoerceToSubfieldMatrix(A, L, K, hKL);
        Aroot, Rroot := Explode(hroot);
        hcomp := [* A0*Aroot, R*Rroot *];
        Append(~comps, [* Qroot, hcomp, incdataroot *]);
    end for;
end for;
return comps;

end intrinsic;
