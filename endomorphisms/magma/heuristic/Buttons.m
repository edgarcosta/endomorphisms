/***
 *  Some useful functions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic HeuristicEndomorphismAlgebra(X::. : Geometric := false, CC := false) -> .
{Returns the endomorphism algebra of X, by default over the base and over QQbar if Geometric is set to true. The first component is the algebra, the second the generators of the endomorphism ring, and the final a string description of the algebra tensored with RR. The second return value is a string description. If CC is set to true, then no algebraization occurs.}

assert ISA(Type(X),Crv) or ISA(Type(X), SECurve);
if CC then
    Geometric := true;
    GeoEndoRep := GeometricEndomorphismRepresentationCC(X);
else
    GeoEndoRep := GeometricEndomorphismRepresentation(X);
end if;

if Geometric then
    EndoAlg, EndoDesc := EndomorphismStructure(GeoEndoRep);
    return EndoAlg, EndoDesc;
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
EndoAlg, EndoDesc := EndomorphismStructure(X`base_endo_rep);
return EndoAlg, EndoDesc;

end intrinsic;


intrinsic HeuristicEndomorphismRing(X::. : Geometric := false, CC := false) -> .
{Returns the endomorphism algebra of X, by default over the base and over QQbar if Geometric is set to true. The first component is the algebra, the second the generators of the endomorphism ring, and the final a string description of the algebra tensored with RR. The second return value is a string description. If CC is set to true, then no algebraization occurs.}

assert ISA(Type(X), Crv) or ISA(Type(X), SECurve);
if CC then
    Geometric := true;
    GeoEndoRep := GeometricEndomorphismRepresentationCC(X);
else
    GeoEndoRep := GeometricEndomorphismRepresentation(X);
end if;

if Geometric then
    EndoAlg, EndoDesc := EndomorphismStructure(GeoEndoRep);
    return Order(Integers(), EndoAlg[2]);
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
EndoAlg, EndoDesc := EndomorphismStructure(X`base_endo_rep);
return Order(Integers(), EndoAlg[2]);

end intrinsic;


intrinsic HeuristicEndomorphismRepresentation(X::. : Geometric := false, CC := false) -> .
{Returns the endomorphism representation of X, by default over the base and over QQbar if Geometric is set to true. If CC is set to true, then no algebraization occurs.}

assert ISA(Type(X),Crv) or ISA(Type(X), SECurve);
if CC then
    Geometric := true;
    GeoEndoRep := GeometricEndomorphismRepresentationCC(X);
else
    GeoEndoRep := GeometricEndomorphismRepresentation(X);
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


intrinsic HeuristicEndomorphismFieldOfDefinition(X::.) -> .
{Returns the field of definition of the endomorphisms of X.}

assert ISA(Type(X), Crv) or ISA(Type(X), SECurve);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
return BaseRing(GeoEndoRep[1][1]);

end intrinsic;


intrinsic HeuristicEndomorphismLattice(X::.) -> .
{Returns the endomorphism lattice of X.}

assert ISA(Type(X), Crv) or ISA(Type(X), SECurve);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
return EndomorphismLattice(GeoEndoRep);

end intrinsic;


intrinsic HeuristicEndomorphismDescription(X::.) -> .
{Returns a description of the endomorphism lattice of X.}

assert ISA(Type(X), Crv) or ISA(Type(X), SECurve);
GeoEndoRep := GeometricEndomorphismRepresentation(X);
return EndomorphismLattice(GeoEndoRep);

end intrinsic;


intrinsic HeuristicIsGL2(X::. : Definition := "Generalized") -> .
{Returns whether or not X is of GL_2-type in the generalized sense (by default) or in the sense of Ribet (if Definition is set to "Ribet").}

assert ISA(Type(X), Crv) or ISA(Type(X), SECurve);
assert Definition in [ "Generalized", "Ribet" ];
if ISA(Type(X), Crv) then
    g := Genus(X);
elif Type(X) eq SECurve then
    g := X`Genus;
end if;
if Definition eq "Generalized" then
    A := HeuristicEndomorphismAlgebra(X)[1];
    if not Dimension(A) eq g then
        return false;
    end if;
    return Center(A) eq A;
elif Definition eq "Ribet" then
    A := HeuristicEndomorphismAlgebra(X)[1];
    if not Dimension(A) eq g then
        return false;
    end if;
    if not Center(A) eq A then
        return false;
    end if;
    return #CentralIdempotents(A) eq 1;
end if;

end intrinsic;


intrinsic HeuristicJacobianFactors(X::. : AllIdems := true, AllPPs := false, ProjToIdem := true, ProjToPP := true) -> .
{Returns factors of the Jacobian of X over the smallest possible fields, together with maps to these factors. Setting AllMaps to true returns multiple entries for a given components in the decomposition together with all possible maps (instead of a single one). Setting AllPPs to true returns multiple entries for a given idempotent, corresponding to the various choices of principal polarization. Setting ProjToIdem to false uses an inclusion instead of a projection when taking idempotents. Setting ProjToPP to false uses an inclusion instead of a projection when making a period matrix principally polarized. If ProjToIdem and ProjToPP are not equal, then right now the algorithm only returns a component, not a corresponding map from or to the Jacobian of X.}

assert ISA(Type(X),Crv) or ISA(Type(X), SECurve);
P := PeriodMatrix(X); gP := #Rows(P);
GeoEndoRep := GeometricEndomorphismRepresentation(X);

/* The upcoming is badly written boilerplate code, but for me it describes the
 * case distinctions (list or lists of lists) fairly */
if not AllIdems and not AllPPs then
    comps := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recs := [ ];
    for comp in comps do
        Q, mor := Explode(comp);
        if #Rows(Q) eq #Rows(P) then
            return [ [* X, [* IdentityMatrix(BaseRing(X), gP), IdentityMatrix(Integers(), 2*gP) *] *] ];
        end if;
        rec := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
        Append(~recs, rec);
    end for;
    return recs;

elif not AllIdems and AllPPs then
    comps := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recss := [ ];
    for comp in comps do
        Q, mor := Explode(comp);
        if #Rows(Q) eq #Rows(P) then
            return [ [ [* X, [* IdentityMatrix(BaseRing(X), gP), IdentityMatrix(Integers(), 2*gP) *] *] ] ];
        end if;
        recs := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
        Append(~recss, recs);
    end for;
    return recss;

elif AllIdems and not AllPPs then
    comptups := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recss := [ ];
    for comptup in comptups do
        recs := [ ];
        for comp in comptup do
            Q, mor := Explode(comp);
            if #Rows(Q) eq #Rows(P) then
                return [ [ [* X, [* IdentityMatrix(BaseRing(X), gP), IdentityMatrix(Integers(), 2*gP) *] *] ] ];
            end if;
            rec := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
            Append(~recs, rec);
        end for;
        Append(~recss, recs);
    end for;
    return recss;

elif AllIdems and AllPPs then
    comptups := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recsss := [ ];
    for comptup in comptups do
        recss := [ ];
        for comp in comptup do
            Q, mor := Explode(comp);
            if #Rows(Q) eq #Rows(P) then
                return [ [ [ [* X, [* IdentityMatrix(BaseRing(X), gP), IdentityMatrix(Integers(), 2*gP) *] *] ] ] ];
            end if;
            recs := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
            Append(~recss, recs);
        end for;
        Append(~recsss, recss);
    end for;
    return recsss;
end if;

end intrinsic;


intrinsic IsogenyInformation(X::. : facinfo := 0) -> .
{Returns homology exponents of isogeny induced by splitting, and tests if it is compatible with the various polarizations. The information in facinfo has to be calculated with AllIdems set to true.}

assert ISA(Type(X),Crv) or ISA(Type(X), SECurve);
if Type(facinfo) eq RngIntElt then
    facinfo := HeuristicJacobianFactors(X);
end if;

Rs := &cat[ [ fac[2][2] : fac in facs ] : facs in facinfo ];
print "OINK";
print facinfo;
print Rs;
gYs := [ #Rows(R) div 2 : R in Rs ];
gX := &+gYs;
EYs := [ StandardSymplecticMatrix(gY) : gY in gYs ];
EX := StandardSymplecticMatrix(gX);

T := VerticalJoin(Rs);
S := SmithForm(ChangeRing(T, Integers()));
D := Diagonal(S);
exps := [ d : d in D ];

EY := DiagonalJoin(EYs);
test := IsRationalMultiple(Transpose(T)*EY*T, EX);
degs := [ Abs((R*EX*Transpose(R))[1, 1 + (#Rows(R) div 2)]) : R in Rs ];

return exps, test, degs;

end intrinsic;


intrinsic CertifiedEndomorphismAlgebra(X::Crv : P0 := 0, Geometric := false, Al := "Cantor", Cheat := false) -> .
{Returns the (certified) endomorphism algebra of X using the base point P0 (if given), by default over the base and over QQbar if Geometric is set to true. The output is the same as for HeuristicEndomorphismAlgebra with the last argument the certificates. Al is either "Cantor" or "Divisor".}

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
print X;

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
