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


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, GalK::List) -> List
{Given a representation of the geometric endomorphism ring and a Galois group
GalK, returns the endomorphism structure over the subfield corresponding to
GalK.}

/* Called Base because it is the version without Sato-Tate */
vprint EndoFind : "";
vprint EndoFind : "Calculating representation over subfield...";
EndoRep := EndomorphismRepresentation(GeoEndoRep, GalK);
vprint EndoFind : "done.";
EndoAlg, EndoDesc := EndomorphismStructure(EndoRep);
EndoData := [* EndoRep, EndoAlg, EndoDesc *];
return EndoData;

end intrinsic;


intrinsic EndomorphismData(GeoEndoRep::SeqEnum, K::Fld, h::Map) -> List
{Given a representation of the geometric endomorphism ring and a field K,
returns the endomorphism structure over K.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return EndomorphismData(GeoEndoRep, GalK);

end intrinsic;


intrinsic EndomorphismDataWithSatoTate(GeoEndoRep::SeqEnum, GalK::List : Shorthand := "") -> List
{Given a representation of the geometric endomorphism ring and Galois group
GalK, returns the endomorphism structure over the subfield corresponding to
GalK. Also calculates Sato-Tate group from Shorthand.}

/* Adds Sato-Tate */
EndoData := EndomorphismData(GeoEndoRep, GalK);
EndoRep, EndoAlg, EndoDesc := Explode(EndoData);
SatoTate := SatoTateGroup(EndoData, GeoEndoRep, GalK : Shorthand := Shorthand);
Append(~EndoAlg, SatoTate); Append(~EndoDesc, SatoTate);
EndoDataWithST := [* EndoRep, EndoAlg, EndoDesc *];
return EndoDataWithST;

end intrinsic;


intrinsic EndomorphismDataWithSatoTate(GeoEndoRep::SeqEnum, K::Fld, h::Map : Shorthand := "") -> List
{Given a representation of the geometric endomorphism ring and a field K,
returns the endomorphism structure over K. Also calculates Sato-Tate group from
Shorthand.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return EndomorphismDataWithSatoTate(GeoEndoRep, GalK : Shorthand := Shorthand);

end intrinsic;


intrinsic EndomorphismStructure(EndoRep::SeqEnum) -> List
{Given a representation EndoRep of an endomorphism ring, returns a description
of the corresponding algebra, ring, and algebra tensored with RR.}

Rs := [ gen[2] : gen in EndoRep ];
vprint EndoFind : "";
vprint EndoFind : "Generators of endomorphism algebra:", Rs;
vprint EndoFind : "Calculating structure...";
/* Creation of relevant algebras */
g := #Rows(Rs[1]) div 2;
/* Ambient matrix algebra, plus generators of the endomorphism ring */
A := Algebra(MatrixRing(Rationals(), 2*g));
GensA := [ A ! Eltseq(R) : R in Rs ];
/* As a subalgebra */
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
/* As an associative algebra */
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

EndoAlg := [* *]; EndoDesc := [* *];
EndoAlgQQ, EndoDescQQ := EndomorphismAlgebraQQ(C);
Append(~EndoAlg, EndoAlgQQ); Append(~EndoDesc, EndoDescQQ);
EndoAlgZZ, EndoDescZZ := EndomorphismAlgebraZZ(C, GensC);
Append(~EndoAlg, EndoAlgZZ); Append(~EndoDesc, EndoDescZZ);
EndoAlgRR, EndoDescRR := EndomorphismAlgebraRR(C, EndoDescQQ);
Append(~EndoAlg, EndoAlgRR); Append(~EndoDesc, EndoDescRR);
vprint EndoFind : "done";
return EndoAlg, EndoDesc;

end intrinsic;


intrinsic EndomorphismAlgebraQQ(C::AlgAss) -> .
{Given an associative algebra C, returns a description of it.}

/* Central decomposition */
Ds := DirectSumDecomposition(C);
EndoDescQQ := [* *];
for D in Ds do
    DescFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := BaseRing(E1);
    E2 := ChangeRing(E1, F);
    FDesc := FieldDescriptionExtra(Polredbestabs(F));
    //FDesc := Eltseq(MinimalPolynomial(F.1));
    //FDesc := [ Integers() ! c : c in FDesc ];

    _, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        if d eq 1 then
            DescFactorQQ := [* "I", FDesc, 1, 1, 1 *];

        elif d eq 2 then
            /* TODO: Depends on genus <= 3:
             * Problem is matrix algebra over division algebra */
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Integers() ! Norm(DQFin);
            if NDQ eq 1 then
                DescFactorQQ := [* "I", FDesc, 1, 1, 2 *];
            elif not IsDefinite(Q) then
                DescFactorQQ := [* "II", FDesc, 2, NDQ, 1 *];
            else
                DescFactorQQ := [* "III", FDesc, 2, NDQ, 1 *];
            end if;

        elif d eq 3 then
            /* TODO: Depends on genus <= 3 */
            DescFactorQQ := [* "I", FDesc, 1, 1, d *];

        else
            DescFactorQQ := [* "I, II or III", FDesc, -1, -1, -1 *];
        end if;

    else
        /* TODO: Depends on genus <= 3 */
        if d le 3 then
            DescFactorQQ := [* "IV", FDesc, 1, 1, d *];
        else
            DescFactorQQ := [* "IV", FDesc, -1, -1, -1 *];
        end if;
    end if;

    Append(~EndoDescQQ, DescFactorQQ);
end for;

return C, EndoDescQQ;

end intrinsic;


intrinsic EndomorphismAlgebraRR(C::AlgAss, EndoDescQQ::List) -> .
{Given an associative algebra C and its description over QQ, returns a
description of the algebra tensored with RR.}
/* TODO: Depends on genus <= 3 */

EndoDescRR := [ ];
for DescFactorQQ in EndoDescQQ do
    AlbertType := DescFactorQQ[1];
    e := #DescFactorQQ[2] - 1;
    d := DescFactorQQ[3];
    m := DescFactorQQ[5];

    if AlbertType eq "I" then
        if m eq 1 then
            str := "RR";
        else
            str := Sprintf("M_%o (RR)", m);
        end if;
        EndoDescRR cat:= [ str : i in [1..e] ];
    elif AlbertType eq "II" then
        EndoDescRR cat:= [ "M_2 (RR)" : i in [1..e] ];
    elif AlbertType eq "III" then
        EndoDescRR cat:= [ "HH" : i in [1..e] ];
    elif AlbertType eq "I, II or III" then
        EndoDescRR cat:= [ "RR, M_2 (RR) or HH" : i in [1..e] ];
    elif AlbertType eq "IV" then
        if m eq 1 then
            str := "CC";
        else
            str := Sprintf("M_%o (CC)", m);
        end if;
        EndoDescRR cat:= [ str : i in [1..(e div 2)] ];
    end if;
end for;
Sort(~EndoDescRR);
return EndoDescRR, EndoDescRR;

end intrinsic;


intrinsic EndomorphismAlgebraZZ(C::AlgAss, GensC::SeqEnum) -> .
{Given an associative algebra C and generators GensC of an order in it, returns
a description of said order.}

/* Calculating index */
OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);

/* Test whether the order is Eichler in a quaternion algebra */
Ds := DirectSumDecomposition(C);
if #Ds eq 1 then
    E1, f1 := AlgebraOverCenter(C);
    //F := ClearDenominator(BaseRing(E1));
    //F := Polredbestabs(F);
    //E2, f2 := ChangeRing(E1, F);
    F := BaseRing(E1);
    E2 := E1;
    test, d := IsSquare(Dimension(E2));

    if d eq 2 and IsQQ(F) then
        test, Q, f3 := IsQuaternionAlgebra(E2);
        if test then
            //f := f1 * f2 * f3;
            f := f1 * f3;
            OO := QuaternionOrder([ f(gen) : gen in GensC ]);
            if IsEichler(OO) then
                return GensC, [ Integers() ! ind, 1 ];
            else
                return GensC, [ Integers() ! ind, 0 ];
            end if;
        end if;
    end if;
end if;
return GensC, [ Integers() ! ind, -1 ];

end intrinsic;


/* The following two functions are not essential */

intrinsic HasGenerator(EndoStruct::List : B := 1) -> BoolElt, .
{Determines whether a single generator for the endomorphism ring exists, and
returns it if it does.}

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

end intrinsic;


intrinsic FewGenerators(EndoStruct::List) -> SeqEnum
{Returns a small generating set of the endomorphism ring corresponding to
EndoStruct.}

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

end intrinsic;
