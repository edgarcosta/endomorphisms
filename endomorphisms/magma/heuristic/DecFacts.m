/***
 *  Decomposing Jacobians
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
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

// TODO: Autoduality (project)
// End result:
//   Give field or subgroup, find isotypical components, roots, and isogenous genus-2 for it,
//   all over the base and the latter over an extension if geometric option is set to true.
//   Join by an overall geometric option that identifies the decomposition field first by using the tangent representation;
//   more detailed info can be obtained by closer inspection, but should not be given by default.


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


intrinsic IsotypicalIdempotents(P::., GeoEndoRep::SeqEnum) -> .
{Returns factors of the Jacobian and a spanning set of idempotents.}

GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep : CalcPic := false);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
C := GeoEndoAlg[1]; idemsC := CentralIdempotents(C);
return [ MatricesFromIdempotent(idemC, GeoEndoData) : idemC in idemsC ];

end intrinsic;


// TODO: We should not algebraize again here and remove the subfield
intrinsic ComponentFromIdempotent(P::., idem::List : CoerceToBase := true, ProjToIdem := true) -> .
{Returns component and map corresponding to idempotent idem. If CoerceToBase is set to true, this is returned over the smallest possible field.}

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Determining component from idempotent, analytic step...";
L := BaseRing(idem[1]);
if ProjToIdem then
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
vprint EndoFind, 2 : "done determining component from idempotent analytically.";

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Determining component from idempotent, algebraic step...";
/* Recalculation to algebraize entries */
test, B := AlgebraizeMatrixExtra(BCC, L);
assert test;
K, hKL := SubfieldExtra(L, Eltseq(B));
incdata := [* L, K, hKL *];

if CoerceToBase then
    B := CoerceToSubfieldMatrix(B, L, K, hKL);
end if;
vprint EndoFind, 2 : "done determining component from idempotent algebraically.";
return Q, [* B, S *], incdata;

end intrinsic;


intrinsic IsotypicalComponents(P::., EndoRep::SeqEnum : CoerceToBase := true, ProjToIdem := true) -> .
{Returns isotypical components Q = B^d of the Jacobian with maps. If CoerceToBase is set to true, these are returned over the smallest possible field.}

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Determining isotypical components...";
idems := IsotypicalIdempotents(P, EndoRep);
comps := [ ];
for idem in idems do
    Q, mor, incdata := ComponentFromIdempotent(P, idem : CoerceToBase := CoerceToBase, ProjToIdem := ProjToIdem);
    Append(~comps, [* Q, mor, incdata *]);
end for;
vprint EndoFind, 2 : "done determining isotypical components.";
return comps;

end intrinsic;


intrinsic SplittingIdempotents(Q::., mor::., incdata::.) -> .
{Returns further idempotents over the smallest field where the isotypical component splits as far as possible.}
// TODO: Now we start over a given base K and recalculate;
//       we find a place in the lattice between K and L where everything shows up for the first time.
//       This is perhaps not yet the smallest field of definition of a single factor!
//       Moreover, all of this should be read off from the lattice description and be made less ad hoc.

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Finding further splitting...";
L, K, hKL := Explode(incdata);

// TODO: Superfluous recalculation of endomorphism algebra
GeoEndoRepCC := GeometricEndomorphismRepresentationCC(Q);
GeoEndoRep := [ ];
vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Algebraizing matrices...";
for tupCC in GeoEndoRepCC do
    test, A := AlgebraizeMatrixExtra(tupCC[1], L); R := tupCC[2];
    Append(~GeoEndoRep, [* A, R *]);
end for;
vprint EndoFind, 3 : "done algebraizing matrices.";

vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Finding geometric endomorphisms...";
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep : CalcPic := false);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
idems_geo := SplittingIdempotentsAlgebra(GeoEndoData);
vprint EndoFind, 3 : "done finding endomorphisms.";

vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Running through lattice...";
/* Find automorphism group (over QQ for now) */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
H := FixedGroupExtra(L, K, hKL);

/* Find subgroups corresponding to fields between L and K */
S := Subgroups(H); Js := [ rec`subgroup : rec in S ];
Js := [ J : J in Js | #J ne 1 ];
if #Js eq 0 then
    vprint EndoFind, 3 : "done running through lattice.";
    vprint EndoFind, 2 : "done finding further splitting.";
    return idems_geo, [* L, L, CanonicalInclusionMap(L, L) *];
end if;

/* Run through lattice and use smallest field where decomposition occurs */
for J in Reverse(Js) do
    gensM := Generators(J); GalM := [* gensM, Gphi *];
    EndoData, hML := EndomorphismData(GeoEndoRep, GalM);
    idems := SplittingIdempotentsAlgebra(EndoData);
    M := BaseRing(idems[1][1]);
    if #idems eq #idems_geo then
        vprint EndoFind, 3 : "done running through lattice.";
        vprint EndoFind, 2 : "done finding further splitting.";
        return idems, [* L, M, hML *];
    end if;
end for;
vprint EndoFind, 3 : "done running through lattice.";
vprint EndoFind, 2 : "done finding further splitting.";
return idems_geo, [* L, L, CanonicalInclusionMap(L, L) *];

end intrinsic;


intrinsic SplittingIdempotentsAlgebra(EndoData::.) -> .
{Returns further decomposition of the algebra defined by EndoData.}
// TODO: Same conditions as in Structure.
//       Moreover, do this per component only.

/* Assert small dimension */
g := #Rows(EndoData[1][1][1]);
assert g le 7;

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
// TODO: Same conditions as in Structure.
//       Moreover, do this per component only.

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
// TODO: Same conditions as in Structure.
//       Moreover, do this per component only.

E1, f1 := AlgebraOverCenter(C);

/* Fields have no idempotents */
if IsCommutative(E1) then
    return true, [ C ! 1 ];
end if;

/* Central algebras of dimension 4 */
// TODO: Unify with next using methods in Structure.m
test_dim, d := IsSquare(Dimension(E1));
if d eq 2 then
    test_quat, Q, f3 := IsQuaternionAlgebra(E1);
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
// TODO: Same conditions as in Structure.
//       Moreover, do this per component only.
/* Catches M3 over QQ or CM, M4 over ZZ or CM, M2 over quaternion algebra */

C := EndoData[2][1]; d := Dimension(C);
dim := 0; goal := 0;
if d in { 9, 18 } then
    dim := 1; goal := 3;
elif d in { 32 } then
    dim := 1; goal := 4;
elif d in { 16 } then
    disc := Abs(Discriminant(MaximalOrder(C)));
    if disc eq 2^32 then
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
    if IsEven(counter) then Bmax +:= 1; else Bmin -:= 1; end if;
    counter +:= 1;
    D := [ Bmin..Bmax ];

    CP := CartesianPower(D, d);
    for tup in CP do
        c := C ! [ x : x in tup ];
        v := Eltseq(c); v2 := Eltseq(c^2);
        test_mult, lambda := IsRationalMultiple(v2, v);
        if test_mult and (lambda ne 0) then
            idem := (1/lambda)*c;
            test_true, idem := IsTrueIdempotent(EndoData, dim, idem);
            if test_true then
                test_seen := MatrixInBasis(idem, idems);
                if not test_seen then
                    Append(~idems, idem);
                end if;
            end if;

            /* Return when enough idempotents were found */
            if #idems eq goal then
                return idems;
            end if;
        end if;
    end for;
end while;

end function;


function IsTrueIdempotent(EndoData, dim, idem)
/* Only accept idempotents whose homology representation has right dimension */

A, R := Explode(MatricesFromIdempotent(idem, EndoData));
if Rank(R) eq 2*dim then
    return true, R;
elif Rank(1 - R) eq 2*dim then
    return true, 1 - R;
end if;
return false;

end function;


intrinsic RootsOfIsotypicalComponent(Q::., mor::., incdata::. : ProjToIdem := true) -> .
{Returns components along with maps over smallest possible field.}
/* We emphatically do not want a single component and multiple maps for it,
 * because that causes us to miss factors that become isogenous only later */

idems, incdataroot := SplittingIdempotents(Q, mor, incdata);

comps := [ ];
for idem in idems do
    Qroot, morroot, _ := ComponentFromIdempotent(Q, idem : CoerceToBase := false, ProjToIdem := ProjToIdem);
    Append(~comps, [* Qroot, morroot, incdataroot *]);
end for;
return comps;

end intrinsic;


// TODO: More logical version: all over given field?
intrinsic SplitComponents(P::., GeoEndoRep::SeqEnum : AllIdems := false, ProjToIdem := true) -> .
{Returns maximal possible splitting of the Jacobian P over the smallest field over which this occurs, plus corresponding projections.}

L := BaseRing(GeoEndoRep[1][1]);
comps := [ ];
vprint EndoFind : "";
vprint EndoFind : "Determining isotypical components...";
// We do not coerce here as this may go too far down
comps_iso := IsotypicalComponents(P, GeoEndoRep : CoerceToBase := false, ProjToIdem := ProjToIdem);
vprint EndoFind : "done determining isotypical components.";
for comp_iso in comps_iso do
    vprint EndoFind : "";
    vprint EndoFind : "Determining isomorphic component...";
    Q, mor, incdata := Explode(comp_iso);
    A, R := Explode(mor);
    comp_roots := RootsOfIsotypicalComponent(Q, mor, incdata : ProjToIdem := ProjToIdem);

    if not AllIdems then
        N := 1;
    else
        N := #comp_roots;
    end if;

    comptup := [ ];
    for i in [1..N] do
        comp_root := comp_roots[i];
        Qroot, morroot, incdataroot := Explode(comp_root);
        L, K, hKL := Explode(incdataroot);
        A0 := CoerceToSubfieldMatrix(A, L, K, hKL);
        Aroot, Rroot := Explode(morroot);
        if ProjToIdem then
            morcomp := [* Aroot*A0, Rroot*R *];
        else
            morcomp := [* A0*Aroot, R*Rroot *];
        end if;

        /* TODO: Remove this sanity check before the Append statement at some point */
        Atest := morcomp[1]; Rtest := morcomp[2];
        ACCtest := EmbedMatrixExtra(Atest); CC := BaseRing(ACCtest);
        if ProjToIdem then
            assert Maximum([ Abs(c) : c in Eltseq(ACCtest*P - Qroot*ChangeRing(Rtest, CC)) ]) le CC`epscomp;
        else
            assert Maximum([ Abs(c) : c in Eltseq(ACCtest*Qroot - P*ChangeRing(Rtest, CC)) ]) le CC`epscomp;
        end if;

        Append(~comptup, [* Qroot, morcomp, incdataroot *]);
    end for;

    if not AllIdems then
        Append(~comps, comptup[1]);
    else
        Append(~comps, comptup);
    end if;
    vprint EndoFind : "";
    vprint EndoFind : "done determining isomorphic component.";
end for;
return comps;

end intrinsic;


// TODO: Some hacks coming up
import "Structure.m": EndomorphismAlgebraQQ;


function ComponentFromIdempotentHack(P, idem);

R := idem[2];
ACC := TangentRepresentation(R, P, P);
Q, mor := ImgIdemp([* ACC, R *], P);
BCC := mor[1]; S := mor[2];
return Q, [* BCC, S *];

end function;


intrinsic DecompositionOverBase(X::.) -> .
{HAX I SAY.}

F := BaseRing(X); P := PeriodMatrix(X); g := #Rows(P);
EndoDesc := HeuristicEndomorphismDescription(X);
EndoRep := X`base_endo_rep;
EndoAlg, EndoDesc := EndomorphismStructure(EndoRep : CalcPic := false);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
C, GensC := Explode(EndoAlg);
EndoAlgQQ, EndoDescQQ, idems := EndomorphismAlgebraQQ(C, GensC, EndoRep : SortResult := false);

// Boring case
if #EndoDescQQ eq 1 then
    e, _, _, _, dim := Explode(EndoDescQQ[1]);
    if e eq 1 then
        return [ [ dim, e ] ], [ ];
    end if;
end if;

facs := [ ]; eqs := [* *];
for i := 1 to #EndoDescQQ do
    // Get factor
    tup := EndoDescQQ[i];
    idem := idems[i]; idemAR := MatricesFromIdempotent(idem, EndoData);
    e, _, _, _, dim := Explode(tup);
    Q, morQ := ComponentFromIdempotentHack(P, idemAR);
    Append(~facs, [ dim, e ]);

    if (dim eq 1) and (e eq 1) then
        if Im(Q[1,1]/Q[1,2]) lt 0 then
            Q := Matrix([ [ Q[1,2], Q[1,1] ] ]);
        end if;
        E := ReconstructCurve(Q, F : Base := true);
        Append(~eqs, E);
        continue;
    end if;

    if (dim eq 2) and (e eq 1) then
        A, R := Explode(morQ);
        EQ := InducedPolarization(StandardSymplecticMatrix(g), R);
        Ts := IsogenousPPLattices(EQ : ProjToPP := true);

        for T in Ts do
            Qp := Q*ChangeRing(Transpose(T), BaseRing(Q));
            Y, _, test := ReconstructCurve(Qp, F : Base := true);
            if test then
                Append(~eqs, Y);
                break;
            end if;
        end for;
        continue;
    end if;

    if (e eq 2) or (e eq 3) then
        GeoEndoRepQ := GeometricEndomorphismRepresentation(Q, F);
        F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRepQ[1][1]));
        EndoRepQ := EndomorphismRepresentation(GeoEndoRepQ, F, h);
        EndoAlgQ, EndoDescQ := EndomorphismStructure(EndoRepQ : CalcPic := false);
        EndoDataQ := [* EndoRepQ, EndoAlgQ, EndoDescQ *];
        idemsQ := SplittingIdempotentsAlgebra(EndoDataQ);

        Qp := ComponentFromIdempotentHack(Q, idemsQ[1]);
        if Im(Qp[1,1]/Qp[1,2]) lt 0 then
            Qp := Matrix([ [ Qp[1,2], Qp[1,1] ] ]);
        end if;
        Ep := ReconstructCurve(Qp, F : Base := true);
        Append(~eqs, Ep);
        continue;
    end if;

    error "Algorithm fall down go boom";
end for;

return facs, eqs;

end intrinsic;


intrinsic IsotypicalField(X::.) -> .
{HAX I SAY.}

EndoRep := GeometricEndomorphismRepresentation(X);
EndoAlg, EndoDesc := EndomorphismStructure(EndoRep : CalcPic := false);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
C, GensC := Explode(EndoAlg);
EndoAlgQQ, EndoDescQQ, idems := EndomorphismAlgebraQQ(C, GensC, EndoRep : SortResult := false);

idemsAR := [ MatricesFromIdempotent(idem, EndoData) : idem in idems ];
L := BaseRing(idemsAR[1][1]);
return SubfieldExtra(L, &cat[ Eltseq(idemAR[1]) : idemAR in idemsAR ]);

end intrinsic;


intrinsic FullDecompositionField(X::.) -> .
{HAX I SAY.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
L := BaseRing(GeoEndoRep[1][1]);
K, h := IsotypicalField(X);
Gp, Gf, Gphi := AutomorphismGroupPari(L);
G := FixedGroupExtra(L, K, h);

Hs := [ rec`subgroup : rec in Subgroups(G) ];
powersums := [ ];
for i in [1..#Hs] do
    H := Hs[i]; gensH := Generators(H);
    Gal := [* gensH, Gphi *];
    EndoData := EndomorphismData(GeoEndoRep, Gal);
    EndoDesc := EndoData[3];
    Append(~powersums, &+[ tup[1] : tup in EndoDesc[2] ]);
end for;

max := Maximum(powersums);
Hs := [ Hs[i] : i in [1..#Hs] | powersums[i] eq max ];
H0 := Hs[#Hs];
test := &and[ (H0 meet H) eq H : H in Hs ];
Hsmax := [ Hmax : Hmax in Hs | &and[ not (((Hmax meet H) eq Hmax) and (Hmax ne H)) : H in Hs ] ];
H0 := &meet(Hsmax);

gensH0 := Generators(H0);
K0, h0 := FixedFieldExtra(L, [ Gphi(gen) : gen in gensH0 ]);
return K0, h0, test;

end intrinsic;


intrinsic DecompositionOverClosure(X::.) -> .
{HAX I SAY.}

P := PeriodMatrix(X); g := #Rows(P);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
L := BaseRing(GeoEndoRep[1][1]);
Gp, Gf, Gphi := AutomorphismGroupPari(L);

K, h := FullDecompositionField(X);
EndoData := EndomorphismData(GeoEndoRep, K, h);
EndoRep, EndoAlg, EndoDesc := Explode(EndoData);
C, GensC := Explode(EndoAlg);
EndoAlgQQ, EndoDescQQ, idems := EndomorphismAlgebraQQ(C, GensC, EndoRep : SortResult := false);

// Boring case
if #EndoDescQQ eq 1 then
    e, _, _, _, dim := Explode(EndoDescQQ[1]);
    if e eq 1 then
        return [ [ dim, e ] ], [ ];
    end if;
end if;

facs := [ ]; eqs := [* *];
for i := 1 to #EndoDescQQ do
    // Get factor
    tup := EndoDescQQ[i];
    idem := idems[i]; idemAR := MatricesFromIdempotent(idem, EndoData);
    e, _, _, _, dim := Explode(tup);
    Q, morQ := ComponentFromIdempotentHack(P, idemAR);
    Append(~facs, [ dim, e ]);

    if (dim eq 1) and (e eq 1) then
        if Im(Q[1,1]/Q[1,2]) lt 0 then
            Q := Matrix([ [ Q[1,2], Q[1,1] ] ]);
        end if;
        E := ReconstructCurve(Q, K : Base := true);
        Append(~eqs, E);
        continue;
    end if;

    if (dim eq 2) and (e eq 1) then
        A, R := Explode(morQ);
        EQ := InducedPolarization(StandardSymplecticMatrix(g), R);
        Ts := IsogenousPPLattices(EQ : ProjToPP := true);

        for T in Ts do
            Qp := Q*ChangeRing(Transpose(T), BaseRing(Q));
            Y, _, test := ReconstructCurve(Qp, K : Base := true);
            if test then
                Append(~eqs, Y);
                break;
            end if;
        end for;
        continue;
    end if;

    if (e eq 2) or (e eq 3) then
        GeoEndoRepQ := GeometricEndomorphismRepresentation(Q, L);
        EndoRepQ := EndomorphismRepresentation(GeoEndoRepQ, K, h);
        EndoAlgQ, EndoDescQ := EndomorphismStructure(EndoRepQ : CalcPic := false);
        EndoDataQ := [* EndoRepQ, EndoAlgQ, EndoDescQ *];
        idemsQ := SplittingIdempotentsAlgebra(EndoDataQ);

        Qp := ComponentFromIdempotentHack(Q, idemsQ[1]);
        if Im(Qp[1,1]/Qp[1,2]) lt 0 then
            Qp := Matrix([ [ Qp[1,2], Qp[1,1] ] ]);
        end if;
        Ep := ReconstructCurve(Qp, K : Base := true);
        Append(~eqs, Ep);
        continue;
    end if;

    error "Algorithm fall down go boom";
end for;

return facs, eqs;

end intrinsic;

