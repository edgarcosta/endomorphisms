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


intrinsic EndomorphismStructureBase(GeoEndoRep::SeqEnum, GalK::List) -> List
{Gives the endomorphism structure over the subfield corresponding to GalK, starting from a list of representations.}

EndoRep := EndomorphismRepresentation(GeoEndoRep, GalK);
EndoAlg, EndoDesc := EndomorphismAlgebraAndDescriptionBase(EndoRep);
EndoStructBase := [* EndoRep, EndoAlg, EndoDesc *];
return EndoStructBase;

end intrinsic;


intrinsic EndomorphismStructureBase(GeoEndoRep::SeqEnum, K::Fld) -> List
{Gives the endomorphism structure over the subfield K, starting from a list of representations.}

L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K);
return EndomorphismStructureBase(GeoEndoRep, GalK);

end intrinsic;


intrinsic EndomorphismStructure(GeoEndoRep::SeqEnum, GalK::List : Shorthand := "") -> List
{Gives the endomorphism structure over the subfield corresponding to GalK, starting from a list of representations.}

EndoStructBase := EndomorphismStructureBase(GeoEndoRep, GalK);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStructBase);
SatoTate := SatoTateGroup(EndoStructBase, GeoEndoRep, GalK : Shorthand := Shorthand);
Append(~EndoAlg, SatoTate); Append(~EndoDesc, SatoTate);
EndoStruct := [* EndoRep, EndoAlg, EndoDesc *];
return EndoStruct;

end intrinsic;


intrinsic EndomorphismStructure(GeoEndoRep::SeqEnum, K::Fld : Shorthand := "") -> List
{Gives the endomorphism structure over the subfield K, starting from a list of representations.}

L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K);
return EndomorphismStructure(GeoEndoRep, GalK);

end intrinsic;


intrinsic EndomorphismAlgebraAndDescriptionBase(EndoRep::SeqEnum) -> List
{Gives the endomorphism structure, starting from a list of representations.}

gensHom := [ gen[2] : gen in EndoRep ];
// Creation of relevant algebras
g := #Rows(gensHom[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring
A := Algebra(MatrixRing(Rationals(), 2*g));
GensA := [ A ! Eltseq(genHom) : genHom in gensHom ];
// As a subalgebra
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
// As an associative algebra
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


intrinsic EndomorphismAlgebraQQBase(C::AlgAss : Optimize := true) -> .
{Decribes the factors of endomorphism algebra.}
// NOTE: Set Optimize to false if this ever becomes a problem (very unlikely)

// Central decomposition
Ds := DirectSumDecomposition(C);
EndoDescQQ := [* *];
for D in Ds do
    DescFactorQQ := [* *];
    E1 := AlgebraOverCenter(D);
    F := BaseRing(E1);
    E2 := ChangeRing(E1, F);
    F := ClearFieldDenominator(F);
    if (Type(F) eq FldNum and Optimize) then
        F := OptimizedRepresentation(F);
        F := ClearFieldDenominator(F);
    end if;
    FDesc := Eltseq(MinimalPolynomial(F.1));
    FDesc := [ Integers() ! c : c in FDesc ];

    test, d := IsSquare(Dimension(E2));
    if IsTotallyReal(F) then
        /* FIXME: This may not be the most logical case distinction */
        if d eq 1 then
            DescFactorQQ := [* "I", FDesc, d, 1 *];

        elif d eq 2 then
            /* FIXME: Depends on genus <= 3 */
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Integers() ! Norm(DQFin);
            if NDQ eq 1 then
                DescFactorQQ := [* "I", FDesc, d, NDQ *];
            elif not IsDefinite(Q) then
                DescFactorQQ := [* "II", FDesc, d, NDQ *];
            else
                DescFactorQQ := [* "III", FDesc, d, NDQ *];
            end if;

        elif d eq 3 then
            /* FIXME: Depends on genus <= 3 */
            DescFactorQQ := [* "I", FDesc, d, 1 *];

        else
            /* FIXME: We do not know what happens otherwise, even when using the
             * extended Albert classification that I have applied. Testing for
             * a matrix ring can be done with
             *     Norm(Discriminant(MaximalOrder(E2)));
             * but that is only a necessary condition; there may be
             * ramification at infinity only, in which case this does not tell
             * enough. For now we get by; this should be addressed with more
             * general functionality for algebras, not by our package. */
            DescFactorQQ := [* "I, II or III", FDesc, d, -1 *];
        end if;

    else
        if d eq 1 then
            DescFactorQQ := [* "IV", FDesc, d, 1 *];
        elif d eq 2 then
            test, Q := IsQuaternionAlgebra(E2);
            DQFin := Discriminant(Q); NDQ := Norm(DQFin);
            DescFactorQQ := [* "IV", FDesc, d, NDQ *];
        elif d eq 3 then
            /* FIXME: Depends on genus <= 3 */
            DescFactorQQ := [* "IV", FDesc, d, 1 *];
        else
            DescFactorQQ := [* "IV", FDesc, d, -1 *];
        end if;
    end if;

    Append(~EndoDescQQ, DescFactorQQ);
end for;

return C, EndoDescQQ;

end intrinsic;


intrinsic EndomorphismAlgebraRRBase(C::AlgAss, EndoDescQQ::List) -> .
{Decribes the factors of endomorphism algebra tensored with RR.}

EndoDescRR := [ ];
for DescFactorQQ in EndoDescQQ do
    AlbertType := DescFactorQQ[1];
    e := #DescFactorQQ[2] - 1;
    d := DescFactorQQ[3];

    if AlbertType eq "I" then
        if d eq 1 then
            str := "RR";
        else
            str := Sprintf("M_%o (RR)", d);
        end if;
        EndoDescRR cat:= [ str : i in [1..e] ];
    elif AlbertType eq "II" then
        EndoDescRR cat:= [ "M_2 (RR)" : i in [1..e] ];
    elif AlbertType eq "III" then
        EndoDescRR cat:= [ "HH" : i in [1..e] ];
    elif AlbertType eq "I, II or III" then
        EndoDescRR cat:= [ "RR, M_2 (RR) or HH" : i in [1..e] ];
    elif AlbertType eq "IV" then
        if d eq 1 then
            str := "CC";
        else
            str := Sprintf("M_%o (CC)", d);
        end if;
        EndoDescRR cat:= [ str : i in [1..(e div 2)] ];
    end if;
end for;
return EndoDescRR, EndoDescRR;

end intrinsic;


intrinsic EndomorphismAlgebraZZBase(C::AlgAss, GensC::SeqEnum : Optimize := true) -> .
{Describes of the endomorphism ring.}

// Calculating index
OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM);

// Test whether Eichler in a quaternion algebra
Ds := DirectSumDecomposition(C);
if #Ds eq 1 then
    E1, f1 := AlgebraOverCenter(C);
    //F := ClearFieldDenominator(BaseRing(E1));
    //if (Type(F) eq FldNum and Optimize) then
    //    F := OptimizedRepresentation(F);
    //    F := ClearFieldDenominator(F);
    //end if;
    //E2, f2 := ChangeRing(E1, F);
    F := BaseRing(E1);
    E2 := E1;
    test, d := IsSquare(Dimension(E2));
    if d eq 2 and Type(F) eq FldRat then
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
