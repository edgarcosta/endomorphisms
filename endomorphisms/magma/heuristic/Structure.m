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


intrinsic EndomorphismStructureBase(GeoEndoRep::SeqEnum, GalK::List, F::Fld) -> List
{Given a representation of the geometric endomorphism ring, a Galois group
GalK, and a base field F, returns the endomorphism structure over the subfield
corresponding to GalK.}

/* Called Base because it is the version without Sato-Tate */
EndoRep := EndomorphismRepresentation(GeoEndoRep, GalK, F);
EndoAlg, EndoDesc := EndomorphismAlgebraAndDescriptionBase(EndoRep);
EndoStructBase := [* EndoRep, EndoAlg, EndoDesc *];
return EndoStructBase;

end intrinsic;


intrinsic EndomorphismStructureBase(GeoEndoRep::SeqEnum, K::Fld, F::Fld) -> List
{Given a representation of the geometric endomorphism ring, a field K, and a
base field F, returns the endomorphism structure over K.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, F);
return EndomorphismStructureBase(GeoEndoRep, GalK, F);

end intrinsic;


intrinsic EndomorphismStructure(GeoEndoRep::SeqEnum, GalK::List, F::Fld : Shorthand := "") -> List
{Given a representation of the geometric endomorphism ring, a Galois group
GalK, and a base field F, returns the endomorphism structure over the subfield
corresponding to GalK. Also calculates Sato-Tate group from Shorthand
representing the geometric endomorphism ring tensored with RR.}

/* Adds Sato-Tate */
EndoStructBase := EndomorphismStructureBase(GeoEndoRep, GalK, F);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStructBase);
SatoTate := SatoTateGroup(EndoStructBase, GeoEndoRep, GalK, F : Shorthand := Shorthand);
Append(~EndoAlg, SatoTate); Append(~EndoDesc, SatoTate);
EndoStruct := [* EndoRep, EndoAlg, EndoDesc *];
return EndoStruct;

end intrinsic;


intrinsic EndomorphismStructure(GeoEndoRep::SeqEnum, K::Fld, F::Fld : Shorthand := "") -> List
{Given a representation of the geometric endomorphism ring, a field K, and a
base field F, returns the endomorphism structure over K. Also calculates
Sato-Tate group from Shorthand representing the geometric endomorphism
ring tensored with RR.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, F);
return EndomorphismStructure(GeoEndoRep, GalK, F : Shorthand := Shorthand);

end intrinsic;


intrinsic EndomorphismAlgebraAndDescriptionBase(EndoRep::SeqEnum) -> List
{Given a representation EndoRep of an endomorphism ring, returns a description
of the corresponding algebra, ring, and algebra tensored with RR.}

gensHom := [ gen[2] : gen in EndoRep ];
/* Creation of relevant algebras */
g := #Rows(gensHom[1]) div 2;
/* Ambient matrix algebra, plus generators of the endomorphism ring */
A := Algebra(MatrixRing(Rationals(), 2*g));
GensA := [ A ! Eltseq(genHom) : genHom in gensHom ];
/* As a subalgebra */
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
/* As an associative algebra */
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

EndoAlg := [* *]; EndoDesc := [* *];
EndoAlgQQ, EndoDescQQ := EndomorphismAlgebraQQBase(C);
Append(~EndoAlg, EndoAlgQQ); Append(~EndoDesc, EndoDescQQ);
EndoAlgZZ, EndoDescZZ := EndomorphismAlgebraZZBase(C, GensC);
Append(~EndoAlg, EndoAlgZZ); Append(~EndoDesc, EndoDescZZ);
EndoAlgRR, EndoDescRR := EndomorphismAlgebraRRBase(C, EndoDescQQ);
Append(~EndoAlg, EndoAlgRR); Append(~EndoDesc, EndoDescRR);
return EndoAlg, EndoDesc;

end intrinsic;


intrinsic EndomorphismAlgebraQQBase(C::AlgAss) -> .
{Given an associative algebra C, returns a description of it.}

/* Central decomposition */
Ds := DirectSumDecomposition(C);
EndoDescQQ := [* *];
for D in Ds do
    DescFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := BaseRing(E1);
    E2 := ChangeRing(E1, F);
    F := ClearFieldDenominator(F);
    if (not IsQQ(F)) then
        F := Polredbestabs(F);
    end if;
    FDesc := Eltseq(MinimalPolynomial(F.1));
    FDesc := [ Integers() ! c : c in FDesc ];

    test, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        /* TODO: This may not be the most logical case distinction */
        if d eq 1 then
            DescFactorQQ := [* "I", FDesc, 1, 1, 1 *];

        elif d eq 2 then
            /* TODO: Depends on genus <= 3 */
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
            /* TODO: We do not know what happens otherwise, even when using the
             * extended Albert classification that I have applied. Testing for
             * a matrix ring can be done with
             *     Norm(Discriminant(MaximalOrder(E2)));
             * but that is only a necessary condition; there may be
             * ramification at infinity only, in which case this does not tell
             * enough. For now we get by; this should be addressed with more
             * general functionality for algebras, not by our package. */
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


intrinsic EndomorphismAlgebraRRBase(C::AlgAss, EndoDescQQ::List) -> .
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
return EndoDescRR, EndoDescRR;

end intrinsic;


intrinsic EndomorphismAlgebraZZBase(C::AlgAss, GensC::SeqEnum) -> .
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
    //F := ClearFieldDenominator(BaseRing(E1));
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
