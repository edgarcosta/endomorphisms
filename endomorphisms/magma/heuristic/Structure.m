/***
 *  Structural description per subfield
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
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


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, GalK::List : CalcPic := true) -> List
{Given a representation of the geometric endomorphism ring and a Galois group
GalK, returns the endomorphism structure over the subfield corresponding to
GalK.}

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Calculating representation over subfield...";
EndoRep, hKL := EndomorphismRepresentation(GeoEndoRep, GalK);
vprint EndoFind, 2 : "done calculating representation over subfield.";
EndoAlg, EndoDesc := EndomorphismStructure(EndoRep : CalcPic := CalcPic);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
return EndoData, hKL;

end intrinsic;


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, K::Fld, h::Map : CalcPic := true) -> List
{Given a representation of the geometric endomorphism ring and a field K,
returns the endomorphism structure over K.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return EndomorphismData(GeoEndoRep, GalK : CalcPic := CalcPic);

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


intrinsic EndomorphismStructure(EndoRep::SeqEnum : CalcPic := true) -> List
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

EndoAlgQQ, EndoDescQQ := EndomorphismAlgebraQQ(C, GensC, EndoRep);
EndoAlgZZ, EndoDescZZ := EndomorphismAlgebraZZ(C, GensC);
EndoDescRR := EndomorphismAlgebraRR(C, EndoDescQQ);
EndoAlg := [* EndoAlgQQ, EndoAlgZZ *];

if CalcPic then
    pic := #RosatiFixedModule(EndoRep);
    EndoDesc := < EndoDescRR, EndoDescQQ, EndoDescZZ, pic >;
else
    EndoDesc := < EndoDescRR, EndoDescQQ, EndoDescZZ >;
end if;
vprint EndoFind, 2 : "done calculating endomorphism structure.";
return EndoAlg, EndoDesc;

end intrinsic;


function EndomorphismAlgebraQQ(C, GensC, EndoRep : SortResult := true)
// Given an associative algebra C, returns a description of it.
// Entry: [ power, dim (D | QQ), field desc, disc (D), dim (factor) ]

// TODO: This uses foul tricks with discriminants of maximal orders.
// Similar discriminantal criteria could be implemented in higher genus.

// The essential case is telling between M_n (B over F) and M_2n (F). It seems
// more worthwhile to find an actual algebra decomposition; the current
// algorithms already cannot tell the discriminant for the division algebra in
// the M_2 (B over F) case in genus 8. Amazingly, everything is OK before that!

// The current version also does not give full determinantal info for QM over a
// proper number field, which happens from genus 4 onward. This could be be
// arranged: One only needs definiteness and the finite primes. However,
// canonically returning descriptions of these is then the next issue. We avoid
// this for now.

// For genus smaller than 8, one can perform this case by hand using the same
// methods as in the algorithms belong and considering the discriminant in
// detail.

/* Central decomposition */
Ds := DirectSumDecomposition(C);
idems := CentralIdempotents(C);

EndoDescQQ := [ ];
for i in [1..#Ds] do
    D := Ds[i]; idem := idems[i];
    E1 := AlgebraOverCenter(D); F := BaseField(E1);
    FDesc := FieldDescriptionExtra(Polredabs(F)); dimF := #FDesc - 1;

    /* Some relative dimensions */
    mdimfac := Rank(MatrixFromIdempotent(C, GensC, idem, EndoRep)) div 2;
    // TODO: This assumption is a bit too strong
    assert mdimfac le 7;
    m2reldimalg := Dimension(E1); test, sqrtm2reldimalg := IsSquare(m2reldimalg);

    DescFactorQQ := < >;
    if IsTotallyReal(F) then
        if sqrtm2reldimalg eq 2 then
            test, Q := IsQuaternionAlgebra(E1);
            if IsMatrixRing(Q) then
                m := 2; disc := 1;
            else
                m := 1; disc := Integers() ! Norm(Discriminant(Q));
                Fac := Factorization(disc); disc := &*[ tup[1] : tup in Fac ];
                if IsOdd(#Fac) then disc *:= -1; end if;
            end if;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        elif sqrtm2reldimalg eq 4 then
            // Here we have the distinction M_4 (QQ) or M_2 (B over QQ)
            discOO := Abs(Discriminant(MaximalOrder(C)));
            if discOO eq 2^32 then
                m := 4; disc := 1;
            else
                m := 2; disc := Abs(discOO) div 2^32;
                Fac := Factorization(disc); disc := &*[ tup[1] : tup in Fac ];
                if IsOdd(#Fac) then disc *:= -1; end if;
            end if;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        elif sqrtm2reldimalg eq 6 then
            // Here we have the distinction M_6 (QQ) or M_3 (B over QQ)
            discOO := Abs(Discriminant(MaximalOrder(C)));
            if discOO eq -2^36 * 3^36 then
                m := 6; disc := 1;
            else
                m := 3; disc := Abs(discOO) div (2^36 * 3^36);
                Fac := Factorization(disc); disc := &*[ tup[1] : tup in Fac ];
                if IsOdd(#Fac) then disc *:= -1; end if;
            end if;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        else
            m := sqrtm2reldimalg; disc := 1;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;
        end if;

    else
        if sqrtm2reldimalg eq 2 then
            test, Q := IsQuaternionAlgebra(E1);
            if IsMatrixRing(Q) then
                m := 2; disc := 1;
            else
                m := 1; disc := Integers() ! Norm(Discriminant(Q));
                Fac := Factorization(disc); disc := &*[ tup[1] : tup in Fac ];
            end if;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;

        else
            m := sqrtm2reldimalg; disc := 1;
            reldimalg := m2reldimalg div m^2; dimfac := mdimfac div m;
            DescFactorQQ := < m, dimF*reldimalg, FDesc, disc, dimfac >;
        end if;
    end if;

    if #DescFactorQQ eq 0 then
        error "All cases in EndomorphismAlgebraQQ fell through!";
        DescFactorQQ := < -1, -1, FDesc, -1, -1 >;
    end if;
    Append(~EndoDescQQ, DescFactorQQ);
end for;
if SortResult then
    EndoDescQQ := Sort(EndoDescQQ);
end if;
return C, EndoDescQQ, idems;

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
/* TODO: Depends on genus (isotypical component) <= 4
         AND additionally genus <= 4 in the CM step */

EndoDescRR := [ ];
for DescFactorQQ in EndoDescQQ do
    if DescFactorQQ[1] eq -1 then
        EndoDescRR cat:= [ <-1, -1> ];
        continue;
    end if;

    m, dimalg, FDesc, disc, dimfac := Explode(DescFactorQQ);
    R := PolynomialRing(Rationals()); F := NumberField(R ! FDesc); e := Degree(F);
    test, d := IsSquare(dimalg div e);
    assert test;

    if IsTotallyReal(F) then
        // Case II/III: Matrix ring over definite/indefinite quaternion algebra
        // over totally real field F
        if d eq 2 then
            if Sign(disc) eq 1 then
                EndoDescRR cat:= [ <d*m, 1> : i in [1..e] ];
                continue;
            else
                EndoDescRR cat:= [ <m, 4> : i in [1..e] ];
                continue;
            end if;

        // Case I: Matrix ring over totally real field F
        else
            EndoDescRR cat:= [ <d*m, 1> : i in [1..e] ];
            continue;
        end if;
    end if;

    // Case IV: Matrix ring over CM field F
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
    test, d := IsSquare(Dimension(E1));

    if d eq 2 and IsQQ(F) then
        test, Q, f3 := IsQuaternionAlgebra(E1);
        if test then
            f := f1 * f3;
            OO := QuaternionOrder([ f(gen) : gen in GensC ]);
            if IsEichler(OO) then
                return GensC, < Integers() ! ind, 1 >;
            else
                return GensC, < Integers() ! ind, 0 >;
            end if;
        end if;
    end if;
end if;
return GensC, < Integers() ! ind, -1 >;

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
