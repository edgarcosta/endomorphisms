/***
 *  Structural description per subfield
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

import "OverField.m": SubgroupGeneratorsUpToConjugacy;
forward EndomorphismAlgebraQQ;
forward MatrixFromIdempotent;
forward EndomorphismAlgebraRR;
forward EndomorphismAlgebraZZ;
forward DecompositionDescription;


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, GalK::List) -> List
{Given a representation of the geometric endomorphism ring and a Galois group
GalK, returns the endomorphism structure over the subfield corresponding to
GalK.}

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Calculating representation over subfield...";
EndoRep, hKL := EndomorphismRepresentation(GeoEndoRep, GalK);
vprint EndoFind, 2 : "done calculating representation over subfield.";
EndoAlg, EndoDesc := EndomorphismStructure(EndoRep);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
return EndoData, hKL;

end intrinsic;


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, K::Fld, h::Map) -> List
{Given a representation of the geometric endomorphism ring and a field K,
returns the endomorphism structure over K.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return EndomorphismData(GeoEndoRep, GalK);

end intrinsic;


//intrinsic EndomorphismDataWithSatoTate(GeoEndoRep::SeqEnum, GalK::List : Shorthand := "") -> List
//{Given a representation of the geometric endomorphism ring and Galois group
//GalK, returns the endomorphism structure over the subfield corresponding to
//GalK. Also calculates Sato-Tate group from Shorthand.}
//
///* Adds Sato-Tate */
//EndoData := EndomorphismData(GeoEndoRep, GalK);
//EndoRep, EndoAlg, EndoDesc := Explode(EndoData);
//SatoTate := SatoTateGroup(EndoData, GeoEndoRep, GalK : Shorthand := Shorthand);
//Append(~EndoAlg, SatoTate); Append(~EndoDesc, SatoTate);
//EndoDataWithST := [* EndoRep, EndoAlg, EndoDesc *];
//return EndoDataWithST;
//
//end intrinsic;


//intrinsic EndomorphismDataWithSatoTate(GeoEndoRep::SeqEnum, K::Fld, h::Map : Shorthand := "") -> List
//{Given a representation of the geometric endomorphism ring and a field K,
//returns the endomorphism structure over K. Also calculates Sato-Tate group from
//Shorthand.}
//
///* Apply previous function after finding a corresponding subgroup */
//L := BaseRing(GeoEndoRep[1][1]);
//GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
//return EndomorphismDataWithSatoTate(GeoEndoRep, GalK : Shorthand := Shorthand);
//
//end intrinsic;


intrinsic EndomorphismStructure(EndoRep::SeqEnum) -> List
{Given a representation EndoRep of an endomorphism ring, return the
corresponding algebra, ring, and algebra tensored with RR, along with
corresponding descriptions.}

Rs := [ gen[2] : gen in EndoRep ]; g := #Rows(Rs[1]) div 2;

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Calculating endomorphism structure...";
/* Creation of relevant algebras */
g := #Rows(Rs[1]) div 2;
/* Ambient matrix algebra, plus generators of the endomorphism ring */
A := Algebra(MatrixRing(Rationals(), 2*g));
GensA := [ A ! Eltseq(R) : R in Rs ];
/* As a subalgebra */
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
/* As an associative algebra */
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

EndoAlg := [* *]; EndoDesc := < >;
EndoAlgQQ, EndoDescQQ := EndomorphismAlgebraQQ(C, GensC, EndoRep);
Append(~EndoAlg, EndoAlgQQ); Append(~EndoDesc, EndoDescQQ);
EndoAlgZZ, EndoDescZZ := EndomorphismAlgebraZZ(C, GensC);
Append(~EndoAlg, EndoAlgZZ); Append(~EndoDesc, EndoDescZZ);
EndoDescRR := EndomorphismAlgebraRR(C, EndoDescQQ); Append(~EndoDesc, EndoDescRR);
pic := #RosatiFixedModule(EndoRep); Append(~EndoDesc, pic);
vprint EndoFind, 2 : "done calculating endomorphism structure.";
return EndoAlg, EndoDesc;

end intrinsic;


function EndomorphismAlgebraQQ(C, GensC, EndoRep)
// Given an associative algebra C, returns a description of it.
/* TODO: Depends on genus <= 3 */
// Entry: [ power, dim (D | QQ), field desc, disc (D), dim (factor) ]

/* Central decomposition */
Ds := DirectSumDecomposition(C);
idems := CentralIdempotents(C);

EndoDescQQ := [ ];
for i in [1..#Ds] do
    D := Ds[i]; idem := idems[i];
    E1 := AlgebraOverCenter(D); F := BaseRing(E1); E2 := ChangeRing(E1, F);
    FDesc := FieldDescriptionExtra(Polredbestabs(F)); dimF := #FDesc - 1;

    /* Some relative dimensions */
    mdimfac := Rank(MatrixFromIdempotent(C, GensC, idem, EndoRep)) div 2;
    m2reldimalg := Dimension(E2);

    if IsTotallyReal(F) then
        if m2reldimalg eq 1 then
            m := 1; disc := 1;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        elif m2reldimalg eq 4 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); disc := Integers() ! Norm(DQFin);
            if disc eq 1 then m := 2; else m := 1; end if;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        elif m2reldimalg eq 9 then
            m := 3; disc := 1;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        else
            DescFactorQQ := < -1, -1, FDesc, -1, -1 >;
        end if;

    else
        if m2reldimalg le 9 then
            test, m := IsSquare(m2reldimalg); disc := 1;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;
        else
            DescFactorQQ := < -1, -1, FDesc, -1, -1 >;
        end if;
    end if;

    Append(~EndoDescQQ, DescFactorQQ);
end for;
return C, Sort(EndoDescQQ);

end function;


function MatrixFromIdempotent(C, GensC, idem, EndoRep)
// Returns homology representation corresponding to idemC

As := [ gen[1] : gen in EndoRep ]; Rs := [ gen[2] : gen in EndoRep ];

idemC := [ Rationals() ! c : c in Eltseq(idem) ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
test, idem := MatrixInBasis(idemC, GensC);
assert test;

idem := Eltseq(idem);
idemR := &+[ idem[i] * Rs[i] : i in [1..#Rs] ];
return idemR;

end function;


function EndomorphismAlgebraRR(C, EndoDescQQ)
// Given an associative algebra C and its description over QQ, returns a
// description of the algebra tensored with RR.

EndoDescRR := [ ];
for DescFactorQQ in EndoDescQQ do
    m, dimalg, FDesc, disc, dimfac := Explode(DescFactorQQ);
    R := PolynomialRing(Rationals()); F := NumberField(R ! FDesc); e := Degree(F);
    reldimalg := dimalg div e;
    reldimdivalg := reldimalg div m^2;

    if DescFactorQQ[1] eq -1 then
        EndoDescRR cat:= [ <-1, -1> : i in [1..e] ];
        continue;
    end if;

    if IsTotallyReal(F) then
        // Case I: Matrix ring over F itself
        if reldimdivalg eq 1 then
            EndoDescRR cat:= [ <m, 1> : i in [1..e] ];
            continue;
        end if;

        // Case II/III: Matrix ring over definite/indefinite quaternion algebra over F
        if reldimdivalg eq 4 then
            if Sign(disc) eq 1 then
                EndoDescRR cat:= [ <2*m, 1> : i in [1..e] ];
                continue;
            else
                EndoDescRR cat:= [ <m, 4> : i in [1..e] ];
                continue;
            end if;
        end if;

        test, d := IsSquare(reldimdivalg);
        EndoDescRR cat:= [ <d*m, 1> : i in [1..e] ];
        continue;
    end if;

    // Case IV: Matrix ring over division algebra over CM field F
    test, d := IsSquare(reldimdivalg);
    EndoDescRR cat:= [ <d*m, 2> : i in [1..(e div 2)] ];
    continue;
end for;
return Sort(EndoDescRR);

end function;


function EndomorphismAlgebraZZ(C, GensC)
// Given an associative algebra C and generators GensC of an order in it,
// returns a description of said order.

/* Calculating index */
OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);

/* Test whether the order is Eichler in a quaternion algebra */
Ds := DirectSumDecomposition(C);
if #Ds eq 1 then
    E1, f1 := AlgebraOverCenter(C);
    F := BaseRing(E1);
    E2 := E1;
    test, d := IsSquare(Dimension(E2));

    if d eq 2 and IsQQ(F) then
        test, Q, f3 := IsQuaternionAlgebra(E2);
        if test then
            f := f1 * f3;
            OO := QuaternionOrder([ f(gen) : gen in GensC ]);
            if IsEichler(OO) then
                return GensC, < Integers() ! ind, 1, DOC >;
            else
                return GensC, < Integers() ! ind, 0, DOC >;
            end if;
        end if;
    end if;
end if;
return GensC, < Integers() ! ind, -1, DOC >;

end function;


/* TODO: The following two functions are of various degrees of redundancy */
function HasGenerator(EndoStruct : B := 1)
// Determines whether a single generator for the endomorphism ring exists, and
// returns it if it does.

g := #Rows(EndoStruct[1][1][1]);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
C, GensC := Explode(EndoAlg);
OC := Order(Integers(), GensC); DOC := Discriminant(OC);
Dom := [-B..B]; CP := CartesianPower(Dom, #GensC);
min := -1;
for tup in CP do
    gen := &+[ tup[i] * GensC[i] : i in [1..#GensC] ];
    minpol := MinimalPolynomial(gen);
    if Degree(minpol) eq #GensC then
        discgen := Abs(Discriminant(minpol));
        //discgen := Abs(Discriminant(EquationOrder(minpol)));
        if discgen eq Abs(DOC) then
            return true, gen;
        else
            if min lt 0 then
                min := discgen;
            end if;
            min := Minimum(min, discgen);
        end if;
    end if;
end for;
return false, min;

end function;


function FewGenerators(EndoStruct)
// Returns a small generating set of the endomorphism ring corresponding to
// EndoStruct.

g := #Rows(EndoStruct[1][1][1]);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
C, GensC := Explode(EndoAlg);
OC := Order(Integers(), GensC); DOC := Discriminant(OC);
Ds := DirectSumDecomposition(C);
if #Ds eq 1 then
    E,f := AlgebraOverCenter(C);
    F := BaseRing(E);
    ZF := Integers(F);
    GensF := [ F ! f(gen) : gen in GensC ];
    for i:=1 to #GensC do
        gens := GensF[1..i];
        if Maximum([ Degree(MinimalPolynomial(gen)) : gen in gens ]) eq #GensC then
            if Abs(Discriminant(sub< ZF | gens >)) eq Abs(DOC) then
                return i;
            end if;
        end if;
    end for;
end if;

end function;
