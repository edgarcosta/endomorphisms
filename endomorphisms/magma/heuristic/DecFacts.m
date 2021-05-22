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
forward SplittingIdempotentsAlgebraStepTwo;
forward FindIdempotentsNaive;
forward IsTrueIdempotent;


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


intrinsic FactorFromIdempotent(P::., idem::List) -> .
{Returns factor and map corresponding to idempotent idem.}

R := idem[2];
ACC := TangentRepresentation(R, P, P);
Q, mor := ImgIdemp([* ACC, R *], P);
BCC := mor[1]; S := mor[2];

L := BaseRing(idem[1]);
test, B := AlgebraizeMatrixExtra(BCC, L);
if not test then error "Failed to algebraize map over base field. Try increasing the precision."; end if;
return Q, [* B, S *];

end intrinsic;


intrinsic IsotypicalIdempotentsOverBase(P::., EndoRep::SeqEnum) -> .
{Returns factors of P and a spanning set of idempotents.}

EndoAlg, EndoDesc := EndomorphismStructure(EndoRep : CalcPic := false);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
C := EndoAlg[1]; idemsC := CentralIdempotents(C);
return [ MatricesFromIdempotent(idemC, EndoData) : idemC in idemsC ];

end intrinsic;


intrinsic IsotypicalFactorsOverBase(P::., EndoRep::SeqEnum) -> .
{Returns isotypical factors Q = B^d of P with maps.}

idems := IsotypicalIdempotentsOverBase(P, EndoRep);
facs := [ ];
for idem in idems do
    Q, mor := FactorFromIdempotent(P, idem);
    Append(~facs, [* Q, mor *]);
end for;
return facs;

end intrinsic;


intrinsic IsotypicalField(X::.) -> .
{Determines the field of definition of the isotypical factors of the Jacobian of X.}

P := PeriodMatrix(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep : CalcPic := false);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];

idemsiso := IsotypicalIdempotentsOverBase(P, GeoEndoRep);
L := BaseRing(GeoEndoRep[1][1]);
Kiso := SubfieldExtra(L, &cat[ Eltseq(idem[1]) : idem in idemsiso ]);
return Kiso;

end intrinsic;


intrinsic SplittingIdempotentsOverBase(P::., EndoRep::SeqEnum) -> .
{Returns an essentially full set of idempotents over the base.}

g := #Rows(P);
if not g le 4 then error "Genus has to be at most 4 for now"; end if;

EndoAlg, EndoDesc := EndomorphismStructure(EndoRep : CalcPic := false);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
return SplittingIdempotentsAlgebra(EndoData);

end intrinsic;


intrinsic SplittingIdempotentsOverClosure(P::., GeoEndoRep::SeqEnum : SmallestField := true) -> .
{Returns an essentially full set of idempotents over the closure.}

g := #Rows(P);
if not g le 4 then error "Genus has to be at most 4 for now"; end if;

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Finding geometric idempotents...";
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep : CalcPic := false);
GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
idemexpsgeo := SplittingIdempotentsAlgebra(GeoEndoData);
vprint EndoFind, 2 : "done finding geometric idempotents.";

L := BaseRing(GeoEndoRep[1][1]);
if not SmallestField then
    vprint EndoFind, 2 : "Return geometric idempotents: SmallestField disabled.";
    return idemexpsgeo, true, CanonicalInclusionMap(L, L);
end if;

if (#idemexpsgeo eq 1) and (idemexpsgeo[1][2] eq 1) then
    vprint EndoFind, 2 : "Return geometric idempotents: Trivial factorization.";
    return idemexpsgeo, true, CanonicalInclusionMap(L, L);
end if;

/* Find automorphism group (over QQ for now) */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
/* Trivial case */
if #Gp eq 1 then
    vprint EndoFind, 2 : "Return geometric idempotents: Trivial automorphism group.";
    return idemexpsgeo, true, CanonicalInclusionMap(L, L);
end if;

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Running through lattice of subgroups..";
Lat := SubgroupLattice(Gp);
Hs := [ Lat[#Lat] ]; H0pairs := [ ];
repeat
    Hsnew := [ ];
    for H in Hs do
        gensH := Generators(H); GalM := [* gensH, Gphi *];
        EndoData, hML := EndomorphismData(GeoEndoRep, GalM : CalcPic := false);
        idemexps := SplittingIdempotentsAlgebra(EndoData);
        if &+[ tup[2] : tup in idemexps ] eq &+[ tupgeo[2] : tupgeo in idemexpsgeo ] then
            Append(~H0pairs, [* H, idemexps, hML *]);
        else
            Hsnew cat:= [ rec`subgroup : rec in MaximalSubgroups(H) ];
        end if;
    end for;
    Hs := [ H : H in Set(Hsnew) ];
until #Hs eq 0;
vprint EndoFind, 2 : "done running through lattice of subgroups.";

max, ind := Maximum([ #tup[1] : tup in H0pairs ]);
Hmax := H0pairs[ind][1]; idemexps := H0pairs[ind][2]; hML := H0pairs[ind][3];
test := IsNormal(Gp, Hmax) and &and[ tup[1] meet Hmax eq tup[1] : tup in H0pairs  ];
return idemexps, test, hML;

end intrinsic;


intrinsic SplittingIdempotentsAlgebra(EndoData::.) -> .
{Returns an essentially full set of idempotents.}

C := EndoData[2][1];
dec := DirectSumDecomposition(C); idemexps := [ ];
for D in dec do
    done := false;
    test1, idemexpsD := SplittingIdempotentsAlgebraStepOne(D, C, EndoData);
    if test1 then
        done := true;
        idemexps cat:= [ [* C ! idemexp[1], idemexp[2] *] : idemexp in idemexpsD ];
    else
        test2, idemspowD := SplittingIdempotentsAlgebraStepTwo(D, C, EndoData);
        if test2 then
            done := true;
            idemexps cat:= [ [* C ! idemexp[1], idemexp[2] *] : idemexp in idemexpsD ];
        end if;
    end if;
    if not done then
        error "All known cases in SplittingIdempotentsAlgebra fell through.";
    end if;
end for;
return [ [* MatricesFromIdempotent(idemexp[1], EndoData), idemexp[2] *] : idemexp in idemexps ];

end intrinsic;


function SplittingIdempotentsAlgebraStepOne(D, C, EndoData)
/* Covers reduced dimension up to 2 */

E1, f1 := AlgebraOverCenter(D);

/* Fields have no idempotents */
if IsCommutative(E1) then
   return true, [ [* D ! 1, 1 *] ];
end if;

/* Central algebras of dimension 4 */
test_dim, d := IsSquare(Dimension(E1));
if d eq 2 then
    test_quat, Q, f3 := IsQuaternionAlgebra(E1);
    test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
    if not test_mat then
        return true, [ [* D ! 1,  1 *] ];
    end if;
    f := f1 * f3 * f4;
    invf := Inverse(f);
    return true, [ [* D ! invf(M ! [1,0,0,0]), 2 *] ];
end if;

return false, [ ];

end function;


function SplittingIdempotentsAlgebraStepTwo(D, C, EndoData);
/* Covers M3 over QQ or CM, M4 over ZZ or CM, M2 over quaternion algebra */

d := Dimension(D);
dim := 0; exp := 0;
if d in { 9, 18 } then
    /* M3 over QQ or CM */
    dim := 1; exp := 3;
elif d in { 32 } then
    /* M4 over CM */
    dim := 1; exp := 4;
elif d in { 16 } then
    /* M4 over ZZ or M2 over quaternion algebra */
    disc := Abs(Discriminant(MaximalOrder(C)));
    if disc eq 2^32 then
        dim := 1; exp := 4;
    else
        dim := 2; exp := 2;
    end if;
end if;

/* Run naive search if applicable */
if dim ne 0 then
    idems := FindIdempotentsNaive(D, C, EndoData, dim, 1);
    return true, [ [* idems[1], exp *] ];
end if;
return false, [ ];

end function;


function FindIdempotentsNaive(D, C, EndoData, dim, exp)

d := Dimension(D);
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
            assert idem^2 eq idem;
            test_true, idem := IsTrueIdempotent(idem, C, EndoData, dim);
            if test_true then
                test_seen := MatrixInBasis(idem, idems);
                if not test_seen then
                    Append(~idems, idem);
                end if;
            end if;

            /* Return when enough idempotents were found */
            if #idems eq exp then
                return idems;
            end if;
        end if;
    end for;
end while;

end function;


function IsTrueIdempotent(idem, C, EndoData, dim)
/* Only accept idempotents whose homology representation has right dimension */

A, R := Explode(MatricesFromIdempotent(C ! idem, EndoData));
if Rank(R) eq 2*dim then
    return true, R;
elif Rank(1 - R) eq 2*dim then
    return true, 1 - R;
end if;
return false;

end function;


intrinsic SplittingFactorsOverBase(P::., EndoRep::SeqEnum) -> .
{Returns an essentially full set of factors over the base.}

idemexps := SplittingIdempotentsOverBase(P, EndoRep);
facs := [ ];
for idemexp in idemexps do
    B, mor := FactorFromIdempotent(P, idemexp[1]);
    Append(~facs, [* B, mor, idemexp[2] *]);
end for;
return facs;

end intrinsic;


intrinsic SplittingFactorsOverClosure(P::., EndoRep::SeqEnum : SmallestField := true) -> .
{Returns an essentially full set of factors over the closure.}

idemexps, test, hML := SplittingIdempotentsOverClosure(P, EndoRep : SmallestField := SmallestField);
facs := [ ];
for idemexp in idemexps do
    B, mor := FactorFromIdempotent(P, idemexp[1]);
    Append(~facs, [* B, mor, idemexp[2] *]);
end for;
return facs, test, hML;

end intrinsic;


intrinsic DecompositionOverBase(P::., EndoRep::SeqEnum : Base := true) -> .
{Returns an essentially full set of factors over the base, along with corresponding curves if they exist. Setting Base to false allows curves over an extension of the base field of EndoRep to describe the factors.}

facs := SplittingFactorsOverBase(P, EndoRep);
if (#facs eq 1) and (facs[1][3] eq 1) then
    return [ [ #Rows(P), 1 ] ], [ ], [ ];
end if;

decbasedesc := [ ]; decbasefacts := [ ]; decbaseeqs := [ ];
for fac in facs do
    B, mor, exp := Explode(fac);
    desc := [ #Rows(B), exp ];
    Append(~decbasedesc, desc);
    fact := [* B, mor *];
    Append(~decbasefacts, fact);
    eqn := ReconstructionFromFactor(P, B, mor : Base := Base);
    Append(~decbaseeqs, eqn);
end for;
return decbasedesc, decbasefacts, decbaseeqs;

end intrinsic;


function ReduceToClasses(S, E)
// Given a list of objects S and an equivalence relation E on S, returns orbits under E.

if #S eq 0 then return []; end if;
if #S eq 1 then return [[S[1]]]; end if;
if #S eq 2 then return E(S[1],S[2]) select [S] else [[S[1]],[S[2]]]; end if;
Ts:=[[S[1]]];
for i:=2 to #S do
    s:=S[i]; sts := true;
    for j:=#Ts to 1 by -1 do // check most recently added entries first in case adjacent objects in S are more likely to be equivalent     (often true)
        if E(s,Ts[j][1]) then Append(~Ts[j],i); sts:=false; break; end if;
    end for;
    if sts then Append(~Ts,[i]); end if;
end for;
return Ts;

end function;


function IsIsogenousSimpleCC(Q1, Q2)

if #Rows(Q1) ne #Rows(Q2) then return false; end if;
GeoHomRep := GeometricHomomorphismRepresentationCC(Q1, Q2);
return #GeoHomRep ne 0;

end function;


intrinsic DecompositionOverClosure(P::., GeoEndoRep::SeqEnum : Base := true, SmallestField := true) -> .
{Returns an essentially full set of factors over the closure, along with corresponding curves if they exist. Setting Base to false allows curves over an extension to describe the factors.}

L := BaseRing(GeoEndoRep[1][1]);
idemsiso := IsotypicalIdempotentsOverBase(P, GeoEndoRep);
Kiso := SubfieldExtra(L, &cat[ Eltseq(idem[1]) : idem in idemsiso ]);

facs, test, hML := SplittingFactorsOverClosure(P, GeoEndoRep : SmallestField := SmallestField);
Kdec := Domain(hML);
Kdecinfo := [* Kdec, test *];

if (#facs eq 1) and (facs[1][3] eq 1) then
    return [ [ #Rows(P), 1 ] ], [ ], [ ], Kiso, Kdecinfo, [ ];
end if;

decgeodesc := [ ]; decgeofacts := [ ]; decgeoeqs := [ ];
for fac in facs do
    B, mor, exp := Explode(fac);
    desc := [ #Rows(B), exp ];
    Append(~decgeodesc, desc);
    fact := [* B, mor *];
    Append(~decgeofacts, fact);
    eqn := ReconstructionFromFactor(P, B, mor : Base := Base);
    Append(~decgeoeqs, eqn);
end for;

if SmallestField then
    merging := [ ];
else
    inds := [1..#decgeofacts];
    function redfunc(i, j); return IsIsogenousSimpleCC(decgeofacts[i][1], decgeofacts[j][1]); end function;
    merging := Sort(ReduceToClasses(inds, redfunc));
end if;

return decgeodesc, decgeofacts, decgeoeqs, Kiso, Kdecinfo, merging;

end intrinsic;


intrinsic DecompositionOverBase(X::. : Base := true) -> .
{Returns root factors B of the Jacobian of X with maps, along with corresponding curves if they exist. Setting Base to false allows curves over an extension of the base field of X to describe the factors.}

P := PeriodMatrix(X);
EndoRep := HeuristicEndomorphismRepresentation(X);
return DecompositionOverBase(P, EndoRep);

end intrinsic;


intrinsic DecompositionOverClosure(X::. : Base := true, SmallestField := false) -> .
{Returns root factors B of the Jacobian of X with maps, along with corresponding curves if they exist. Setting Base to false allows curves over an extension of the endomorphism field of X to describe the factors.}

P := PeriodMatrix(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
return DecompositionOverClosure(P, GeoEndoRep : SmallestField := SmallestField);

end intrinsic;
