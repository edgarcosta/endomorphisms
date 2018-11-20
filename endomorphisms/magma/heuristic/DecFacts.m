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
forward SplittingIdempotentsAlgebraStepOne;
forward SplittingIdempotentsAlgebraStepOneSubstep;
forward SplittingIdempotentsAlgebraStepTwo;
forward FindIdempotentsStupid;
forward IsTrueIdempotent;

/* TODO: Too much information is recalculated: this should be reproved by true abelian functionality if possible */


intrinsic IsotypicalIdempotents(P::., GeoEndoRep::SeqEnum) -> .
{Returns factors of the Jacobian and a spanning set of idempotents.}

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


intrinsic ComponentFromIdempotent(P::., idem::List : CoerceToBase := true, ProjOrInc := "Proj") -> .
{Returns component and map corresponding to idempotent idem. If CoerceToBase is set to true, this is returned over the smallest possible field.}

vprint EndoFind : "";
vprint EndoFind : "Determining component from idempotent, analytic step...";
L := BaseRing(idem[1]);
if ProjOrInc eq "Proj" then
    R := idem[2];
    ACC := TangentRepresentation(R, P, P);
    //Q, mor := ImgProj([* ACC, R *], P, P);
    /* TODO: This line is needed because of a segfault issue in later Magma versions */
    Q, mor := ImgIdemp([* ACC, R *], P);
else
    R := 1 - idem[2];
    ACC := TangentRepresentation(R, P, P);
    Q, mor := Ker0([* ACC, R *], P, P);
end if;
BCC := mor[1]; S := mor[2];
vprint EndoFind : "done determining component from idempotent analytically.";

vprint EndoFind : "";
vprint EndoFind : "Determining component from idempotent, algebraic step...";
/* Recalculation to algebraize entries */
test, B := AlgebraizeMatrix(BCC, L);
assert test;
K, hKL := SubfieldExtra(L, Eltseq(B));
incdata := [* L, K, hKL *];
if CoerceToBase then
    B := CoerceToSubfieldMatrix(B, L, K, hKL);
end if;
vprint EndoFind, 2 : "";
vprint EndoFind : "done determining component from idempotent algebraically.";
return Q, [* B, S *], incdata;

end intrinsic;


intrinsic IsotypicalComponents(P::., EndoRep::SeqEnum : CoerceToBase := true, ProjOrInc := "Proj") -> .
{Returns isotypical components Q = B^d of the Jacobian with maps. If CoerceToBase is set to true, these are returned over the smallest possible field.}

vprint EndoFind : "";
vprint EndoFind : "Determining isotypical components...";
idems := IsotypicalIdempotents(P, EndoRep);
comps := [ ];
for idem in idems do
    Q, mor, incdata := ComponentFromIdempotent(P, idem : CoerceToBase := CoerceToBase, ProjOrInc := ProjOrInc);
    Append(~comps, [* Q, mor, incdata *]);
end for;
vprint EndoFind, 2 : "";
vprint EndoFind : "done determining isotypical components.";
return comps;

end intrinsic;


intrinsic SplittingIdempotents(Q::., mor::., incdata::.) -> .
{Returns further idempotents over the smallest field where the isotypical component splits as far as possible.}

L, K, hKL := Explode(incdata);
/* Recalculate endomorphism algebra over known field (as mentioned above, this is stupid) */
GeoEndoRepCC := GeometricEndomorphismRepresentationCC(Q);
GeoEndoRep := [ ];
vprint EndoFind : "";
vprint EndoFind : "Algebraizing matrices...";
for tupCC in GeoEndoRepCC do
    test, A := AlgebraizeMatrix(tupCC[1], L); R := tupCC[2];
    Append(~GeoEndoRep, [* A, R *]);
end for;
vprint EndoFind, 2 : "";
vprint EndoFind : "done algebraizing matrices.";
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
idems_geo := SplittingIdempotentsAlgebra(GeoEndoData);

/* Find automorphism group (over QQ for now) */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
H := FixedGroupExtra(L, K, hKL);

/* Find subgroups corresponding to fields between L and K */
S := Subgroups(H); Js := [ rec`subgroup : rec in S ];
Js := [ J : J in Js | #J ne 1 ];
if #Js eq 0 then
    return idems_geo, [* L, L, CanonicalInclusionMap(L, L) *];
end if;

/* Run through lattice and use smallest field where decomposition occurs */
for J in Reverse(Js) do
    gensM := Generators(J); GalM := [* gensM, Gphi *];
    EndoData, hML := EndomorphismData(GeoEndoRep, GalM);
    idems := SplittingIdempotentsAlgebra(EndoData);
    M := BaseRing(idems[1][1]);
    if #idems eq #idems_geo then
        return idems, [* L, M, hML *];
    end if;
end for;
return idems_geo, [* L, L, CanonicalInclusionMap(L, L) *];

end intrinsic;


intrinsic SplittingIdempotentsAlgebra(EndoData::.) -> .
{Returns further decomposition of the algebra defined by EndoData.}

/* Assert small dimension */
g := #Rows(EndoData[1][1][1]);
assert g le 4;

/* Try successively more complicated methods */
test1, idems := SplittingIdempotentsAlgebraStepOne(EndoData);
if test1 then
    return [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
end if;
test2, idems := SplittingIdempotentsAlgebraStepTwo(EndoData);
if test2 then
    return [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
end if;
print "All known cases in SplittingIdempotentsAlgebra fell through, returning unit element.";
return [ MatricesFromIdempotent(EndoData[2][1] ! 1, EndoData) : idem in idems ];

end intrinsic;


function SplittingIdempotentsAlgebraStepOne(EndoData)

C := EndoData[2][1];
idems := [ ];
for D in DirectSumDecomposition(C) do
    test, idemsD := SplittingIdempotentsAlgebraStepOneSubstep(D);
    if not test then
        return false, [ ];
    end if;
    idems cat:= [ C ! idem : idem in idemsD ];
end for;
return true, idems;

end function;


function SplittingIdempotentsAlgebraStepOneSubstep(C)

E1, f1 := AlgebraOverCenter(C);
/* TODO: Optionally, we could improve and change ring here, but that has
 * undocumented behavior and does not accept field morphisms, so not for now */
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
    f := f1 * f3 * f4;
    invf := Inverse(f);
    return true, [ C ! invf(M ! [1,0,0,0]), C ! invf(M ! [0,0,0,1]) ];
end if;
return false, [ ];

end function;


function SplittingIdempotentsAlgebraStepTwo(EndoData);
/* Catches M3 over QQ or CM, M4 over ZZ or CM, M2 over quaternion algebra */

C := EndoData[2][1]; d := Dimension(C);
dim := 0; goal := 0;
if d in { 9, 18 } then
    dim := 1; goal := 3;
elif d in { 32 } then
    dim := 1; goal := 4;
elif d in { 16 } then
    disc := Abs(Discriminant(MaximalOrder(C)));
    if disc eq 4294967296 then
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
        c := C ! [ x : x in tup ];
        v := Eltseq(c); v2 := Eltseq(c^2);
        test_mult, lambda := IsMultiple(v2, v);
        if test_mult and (lambda ne 0) then
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
if Rank(R) eq 2*dim then
    return true;
end if;
return false;

end function;


intrinsic RootsOfIsotypicalComponent(Q::., mor::., incdata::. : ProjOrInc := "Proj") -> .
{Returns components along with maps over smallest possible field.}
/* We emphatically do not want a single component and multiple maps for it,
 * because that causes us to miss factors that become isogenous only later */

idems, incdataroot := SplittingIdempotents(Q, mor, incdata);

comps := [ ];
for idem in idems do
    /* No need to coerce since we are already over smallest possible field; in
     * fact the field might even drop too much */
    Qroot, morroot, _ := ComponentFromIdempotent(Q, idem : CoerceToBase := false, ProjOrInc := ProjOrInc);
    Append(~comps, [* Qroot, morroot, incdataroot *]);
end for;
return comps;

end intrinsic;


intrinsic SplitComponents(P::., GeoEndoRep::SeqEnum : ProjOrInc := "Proj", AllMaps := false) -> .
{Returns maximal possible splitting of the Jacobian P over the smallest field over which this occurs, plus corresponding projections.}

L := BaseRing(GeoEndoRep[1][1]);
comps := [ ];
comps_iso := IsotypicalComponents(P, GeoEndoRep : CoerceToBase := false, ProjOrInc := ProjOrInc);
for comp_iso in comps_iso do
    Q, mor, incdata := Explode(comp_iso);
    A, R := Explode(mor);
    comp_roots := RootsOfIsotypicalComponent(Q, mor, incdata : ProjOrInc := ProjOrInc);
    if AllMaps then
        N := #comp_roots;
    else
        N := 1;
    end if;
    for i in [1..N] do
        comp_root := comp_roots[i];
        Qroot, morroot, incdataroot := Explode(comp_root);
        L, K, hKL := Explode(incdataroot);
        A0 := CoerceToSubfieldMatrix(A, L, K, hKL);
        Aroot, Rroot := Explode(morroot);
        if ProjOrInc eq "Proj" then
            morcomp := [* Aroot*A0, Rroot*R *];
        else
            morcomp := [* A0*Aroot, R*Rroot *];
        end if;

        /* TODO: Remove this sanity check before the Append statement at some point */
        Atest := morcomp[1]; Rtest := morcomp[2];
        ACCtest := EmbedMatrixExtra(Atest, K`iota); CC := BaseRing(ACCtest);
        assert Maximum([ Abs(c) : c in Eltseq(ACCtest*P - Qroot*ChangeRing(Rtest, CC)) ]) le CC`epscomp;

        Append(~comps, [* Qroot, morcomp, incdataroot *]);
    end for;
end for;
return comps;

end intrinsic;
