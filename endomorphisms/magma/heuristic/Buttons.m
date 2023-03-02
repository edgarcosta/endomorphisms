/***
 *  Some useful functions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic CurveExtra(X::Crv : prec := 100) -> .
{This converts a curve to be over a NumberFieldExtra if it is not given in that way.}

assert ISA(Type(X), Crv);
F := BaseRing(X);
if Type(F) eq RngInt then
    X := ChangeRing(X, Rationals());
    F := BaseRing(X);
end if;
if assigned F`base then
    return X;
end if;

if IsQQ(F) then
    K := RationalsExtra(prec);
    h := hom< F -> K | >;
    return ChangeRingCurve(X, h);
end if;

K := NumberFieldExtra(DefiningPolynomial(F) : prec := prec);
test, h := IsIsomorphic(F, K);
assert test;
return ChangeRingCurve(X, h);

end intrinsic;


intrinsic HeuristicEndomorphismAlgebra(X::. : Geometric := false, CC := false, UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns the abstract endomorphism algebra of the Jacobian of X, by default over the base and over QQbar if Geometric is set to true. If CC is set to true, the algebra over QQbar is determined without performing any further algebraization. Notational conventions are specified in the documentation folder.}

assert ISA(Type(X), Crv);
X := CurveExtra(X);
if CC then
    Geometric := true;
    GeoEndoRep := GeometricEndomorphismRepresentationCC(X);
else
    GeoEndoRep := GeometricEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
end if;

if Geometric then
    EndoAlg, EndoDesc := EndomorphismStructure(GeoEndoRep);
    return EndoDesc, EndoAlg[1], EndoAlg[2];
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
EndoAlg, EndoDesc := EndomorphismStructure(X`base_endo_rep);
/* Description, algebra, ring */
return EndoDesc, EndoAlg[1], EndoAlg[2];

end intrinsic;


function HumanizeAlt(desc);
/* Make stuff readable */

R<t> := PolynomialRing(Rationals());

str := "";
alg := desc[2];
str cat:= "Endomorphism algebra (one entry per isotypical factor):\n";
for fac in alg do
    m, dimD, pol, discD, dimA := Explode(fac);
    dimK := #pol - 1;
    assert (dimD mod dimK) eq 0;
    indDsq := dimD div dimK;
    test, indD := IsSquare(indDsq);
    assert test;

    if m gt 1 then
        str cat:= (Sprint(m) cat "th power of ");
    end if;
    str cat:= ("abelian variety of dimension " cat Sprint(dimA)) cat " ";
    str cat:= "with endomorphism algebra ";
    if indD ne 1 then
        str cat:= ("a division algebra of index " cat Sprint(indD) cat " and discriminant " cat Sprint(discD) cat " over ");
    end if;
    if #pol eq 2 then
        str cat:= "QQ\n";
    else
        str cat:= ("the number field defined by " cat Sprint(R ! pol)) cat "\n";
    end if;
end for;

ind, test := Explode(desc[3]);
str cat:= "\n";
str cat:= "Endomorphism ring:\n";
if test eq -1 then
    str cat:= "A subring of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra\n";
elif test eq 0 then
    str cat:= "A non-Eichler order of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra\n";
elif test eq 1 then
    str cat:= "An Eichler order of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra\n";
end if;

algRR := desc[1];
str cat:= "\n";
str cat:= "Endomorphism algebra tensored with RR:\n";
descs := [ "RR", "CC", "undef", "HH" ];
for i := 1 to #algRR do
    fac := algRR[i];
    m, desc := Explode(fac);
    if m eq 1 then
        str cat:= Sprint(descs[desc]);
    else
        str cat:= "M_" cat Sprint(m) cat " (" cat Sprint(descs[desc]) cat ")";
    end if;
    if i ne #algRR then
        str cat:= " x ";
    end if;
end for;

return str;

end function;


function Humanize(desc);
/* Make stuff readable */

R<t> := PolynomialRing(Rationals());

str := "";
alg := desc[2];
str cat:= "Endomorphism algebra:\n";
for i := 1 to #alg do
    fac := alg[i];
    m, dimD, pol, discD, dimA := Explode(fac);
    dimK := #pol - 1;
    assert (dimD mod dimK) eq 0;
    indDsq := dimD div dimK;
    test, indD := IsSquare(indDsq);
    assert test;

    if m gt 1 then
        str cat:= "M_" cat Sprint(m) cat " (";
    end if;
    if indD ne 1 then
        assert indD eq 2;
        str cat:= "B_" cat Sprint(discD) cat " (";
    end if;
    if #pol eq 2 then
        str cat:= "QQ";
    else
        str cat:= "QQ [t] / (" cat Sprint(R ! pol) cat ")";
    end if;
    if indD ne 1 then
        str cat:= ")";
    end if;
    if m gt 1 then
        str cat:= ")";
    end if;

    if i ne #alg then
        str cat:= " x ";
    end if;
end for;

ind, test := Explode(desc[3]);
str cat:= "\n";
str cat:= "Endomorphism ring:\n";
if test eq -1 then
    str cat:= "A subring of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra";
elif test eq 0 then
    str cat:= "A non-Eichler order of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra";
elif test eq 1 then
    str cat:= "An Eichler order of index " cat Sprint(ind) cat " in the maximal order of the endomorphism algebra";
end if;

algRR := desc[1];
str cat:= "\n";
str cat:= "Endomorphism algebra tensored with RR:\n";
descs := [ "RR", "CC", "undef", "HH" ];
for i := 1 to #algRR do
    fac := algRR[i];
    m, desc := Explode(fac);
    if m eq 1 then
        str cat:= Sprint(descs[desc]);
    else
        str cat:= "M_" cat Sprint(m) cat " (" cat Sprint(descs[desc]) cat ")";
    end if;
    if i ne #algRR then
        str cat:= " x ";
    end if;
end for;

return str;

end function;


intrinsic HeuristicEndomorphismDescription(X::. : Geometric := false, CC := false, UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns a readable description of the endomorphism algebra of the Jacobian of X, by default over the base and over QQbar if Geometric is set to true. If CC is set to true, the algebra over QQbar is determined without performing any further algebraization.}

return Humanize(HeuristicEndomorphismAlgebra(X : Geometric := Geometric, CC := CC, UpperBound := UpperBound, DegreeDivides := DegreeDivides));

end intrinsic;


intrinsic HeuristicEndomorphismRepresentation(X::. : Geometric := false, CC := false, UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns the endomorphism representation of X the Jacobian of X, by default over the base and over QQbar if Geometric is set to true. If CC is set to true, then no algebraization occurs. Notational conventions are specified in the documentation folder.}

assert ISA(Type(X),Crv);
X := CurveExtra(X);
if CC then
    Geometric := true;
    GeoEndoRep := GeometricEndomorphismRepresentationCC(X);
else
    GeoEndoRep := GeometricEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
end if;

if Geometric then
    return GeoEndoRep;
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
return X`base_endo_rep;

end intrinsic;


intrinsic HeuristicEndomorphismFieldOfDefinition(X::. : UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns the field of definition of the endomorphism algebra of the Jacobian of X.}

assert ISA(Type(X), Crv);
X := CurveExtra(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
return BaseRing(GeoEndoRep[1][1]);

end intrinsic;


intrinsic HeuristicEndomorphismLattice(X::. : UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns an encoded description of the endomorphism lattice of the Jacobian of X.}

assert ISA(Type(X), Crv);
X := CurveExtra(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
return EndomorphismLattice(GeoEndoRep);

end intrinsic;


intrinsic HeuristicEndomorphismLatticeDescription(X::. : UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns a description of the endomorphism lattice of the Jacobian of X.}

assert ISA(Type(X), Crv);
X := CurveExtra(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
Lat := EndomorphismLattice(GeoEndoRep)[3];
R<t> := PolynomialRing(Rationals());

str := "";
for i := 1 to #Lat do
    entry := Lat[i];
    str cat:= "========\n";
    if entry[2] eq [-1, 1] then
        str cat:= "Over QQ:\n";
    else
        str cat:= "Over QQ [t] / (" cat Sprint(R ! entry[2]) cat "):\n";
    end if;
    str cat:= "--------\n";
    str cat:= Humanize(< entry[3], entry[4], entry[5] >);
    if i ne #Lat then
        str cat:= "\n";
    end if;
end for;
return str;

end intrinsic;


intrinsic HeuristicIsGL2(X::. : Geometric := false, Definition := "Exact", UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns whether or not the Jacobian of X is of GL_2-type in the exact sense (by default) or in the general sense (if Definition is set to "General"; this option is specific to genus 2 for the moment). By default, this question is posed over the base ring, but setting "Geometric" to "true" switches to the algebraic closure}

assert ISA(Type(X), Crv);
assert Definition in [ "Exact", "General" ];
X := CurveExtra(X); g := Genus(X);
if Definition eq "Exact" then
    desc, A := HeuristicEndomorphismAlgebra(X : Geometric := Geometric, UpperBound := UpperBound, DegreeDivides := DegreeDivides);
    if not Dimension(A) eq g then
        return false;
    end if;
    return Center(A) eq A;
elif Definition eq "General" then
    assert g eq 2;
    desc, A := HeuristicEndomorphismAlgebra(X : Geometric := Geometric, UpperBound := UpperBound, DegreeDivides := DegreeDivides);
    return Dimension(A) mod g eq 0;
end if;

end intrinsic;


intrinsic HeuristicDecomposition(X::. : UpperBound := 16, DegreeDivides := Infinity(), SmallestField := true) -> .
{Returns a description of the decomposition of the Jacobian of X. The first entry describes the isotypical field, the second entry describes the decomposition field and whether this field is minimal, the third entry describes decomposition over the base (dimensions and powers followed by equations), and the fourth entry describes the decomposition over the decomposition field (dimensions and powers followed by equations).}

assert ISA(Type(X),Crv);
X := CurveExtra(X);
/* The next line yields a precomputation of the endomorphisms algebra */
GeoEndoRep := HeuristicEndomorphismRepresentation(X : Geometric := true, UpperBound := UpperBound, DegreeDivides := DegreeDivides);
EndoRep := HeuristicEndomorphismRepresentation(X : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
if BaseRing(GeoEndoRep[1][1]) eq BaseRing(EndoRep[1][1]) then
    decgeodesc, decgeofacts, decgeoeqs, Kiso, Kdecinfo, merging := DecompositionOverClosure(X);
    return [* Kiso, Kdecinfo, [* decgeodesc, decgeoeqs *], [* decgeodesc, decgeoeqs, merging *] *];
end if;
decbasedesc, decbasefacts, decbaseeqs := DecompositionOverBase(X);
decgeodesc, decgeofacts, decgeoeqs, Kiso, Kdecinfo, merging := DecompositionOverClosure(X : SmallestField := SmallestField);
return [* Kiso, Kdecinfo, [* decbasedesc, decbaseeqs *], [* decgeodesc, decgeoeqs, merging *] *];

end intrinsic;


intrinsic HeuristicDecompositionFactors(X::. : UpperBound := 16, DegreeDivides := Infinity(), SmallestField := true) -> .
{Returns the factors of the Jacobian of X over the smallest field where they all appear.}

return [ tup[1] : tup in HeuristicDecomposition(X)[3][2] ];

end intrinsic;


intrinsic CertifiedEndomorphismAlgebra(X::Crv : P0 := 0, Geometric := false, Al := "Cantor", Cheat := false) -> .
{Returns the (certified) endomorphism algebra of X using the base point P0 (if given), by default over the base and over QQbar if Geometric is set to true. The output is the same as for HeuristicEndomorphismAlgebra with the last argument the certificates. Al is either "Cantor" or "Divisor".}

assert ISA(Type(X),Crv);
X := CurveExtra(X);
F := BaseRing(X);
if IsHyperelliptic(X) and Cheat then
    f, h := HyperellipticPolynomials(X);
    g := 4*f + h^2;
    g /:= LeadingCoefficient(g);
    if IsEven(Degree(g)) then
        X := HyperellipticCurve(g);
    else
        R<x> := Parent(g);
        e := Degree(g) div 2;
        n := 0;
        repeat
            geven := Numerator(Evaluate(g, 1/(x - n)));
            geven /:= LeadingCoefficient(geven);
            n +:= 1;
        until Degree(geven) mod 2 eq 0;
        X := HyperellipticCurve(geven);
    end if;
end if;
vprint EndoFind : "";
vprint EndoFind : "Curve after transformation to standard form:";
vprint EndoFind : X;

GeoEndoRepCC := HeuristicEndomorphismRepresentation(X : Geometric := true, CC := true);
if not VerifySaturated(GeoEndoRepCC, X`period_matrix) then
    return false, "Not saturated";
end if;
EndoRep := HeuristicEndomorphismRepresentation(X : Geometric := Geometric);

fss := [* *];
for rep in EndoRep do
    test, fs := Correspondence(X, X, rep : P := P0, Q := P0, Al := Al);
    if not test then
        return false, rep;
    end if;
    Append(~fss, rep cat [* fs *]);
end for;
return true, fss;

end intrinsic;
