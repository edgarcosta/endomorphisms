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


intrinsic HeuristicEndomorphismFieldOfDefinition(X::Crv) -> .
{Returns the field of definition of the endomorphisms of X.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
return BaseRing(GeoEndoRep[1][1]);

end intrinsic;


intrinsic HeuristicEndomorphismAlgebra(X::Crv : Geometric := false) -> .
{Returns the endomorphism algebra of X. The first component is the algebra, the second the generators of the endomorphism ring, and the final a string description of the algebra tensored with RR.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
if Geometric then
    EndoAlg, EndoDesc := EndomorphismStructure(GeoEndoRep);
    return EndoAlg;
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
EndoAlg, EndoDesc := EndomorphismStructure(X`base_endo_rep);
return EndoAlg;

end intrinsic;


intrinsic HeuristicEndomorphismAlgebraDescription(X::Crv : Geometric := false) -> .
{Returns a string description of the endomorphism algebra of X.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
if Geometric then
    EndoAlg, EndoDesc := EndomorphismStructure(GeoEndoRep);
    return EndoDesc;
end if;
if not assigned X`base_endo_rep then
    F, h := InclusionOfBaseExtra(BaseRing(GeoEndoRep[1][1]));
    X`base_endo_rep := EndomorphismRepresentation(GeoEndoRep, F, h);
end if;
EndoAlg, EndoDesc := EndomorphismStructure(X`base_endo_rep);
return EndoDesc;

end intrinsic;


intrinsic HeuristicEndomorphismLattice(X::Crv) -> .
{Returns the endomorphism lattice of X.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
return EndomorphismLattice(GeoEndoRep);

end intrinsic;


intrinsic HeuristicIsGL2Ribet(X::Crv) -> .
{Returns whether or not X is of GL_2-type in the sense of Ribet.}

A := HeuristicEndomorphismAlgebra(X)[1];
if not Dimension(A) eq Genus(X) then
    return false;
end if;
if not Center(A) eq A then
    return false;
end if;
return #CentralIdempotents(A) eq 1;

end intrinsic;


intrinsic HeuristicIsGL2Generalized(X::Crv) -> .
{Returns whether or not X is of GL_2-type in the generalized sense of Booker--Sijsling--Sutherland--Voight--Yasaki.}

A := HeuristicEndomorphismAlgebra(X)[1];
if not Dimension(A) eq Genus(X) then
    return false;
end if;
return Center(A) eq A;

end intrinsic;


intrinsic HeuristicJacobianFactors(X::Crv : AllIdems := false, AllPPs := false, ProjToIdem := true, ProjToPP := true) -> .
{Returns factors of the Jacobian of X over the smallest possible fields, together with maps to these factors. Setting AllMaps to true returns multiple entries for a given components in the decomposition together with all possible maps (instead of a single one). Setting AllPPs to true returns multiple entries for a given idempotent, corresponding to the various choices of principal polarization. Setting ProjToIdem to false uses an inclusion instead of a projection when taking idempotents. Setting ProjToPP to false uses an inclusion instead of a projection when making a period matrix principally polarized. If ProjToIdem and ProjToPP are not equal, then right now the algorithm only returns a component, not a corresponding map from or to the Jacobian of X.}

P := PeriodMatrix(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X);

/* The upcoming is badly written boilerplate code, but for me it describes the
 * case distinctions (list or lists of lists) fairly */
if not AllIdems and not AllPPs then
    comps := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recs := [ ];
    for comp in comps do
        Q, mor := Explode(comp);
        rec := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
        Append(~recs, rec);
    end for;
    return recs;

elif not AllIdems and AllPPs then
    comps := SplitComponents(P, GeoEndoRep : AllIdems := AllIdems, ProjToIdem := ProjToIdem);
    recss := [ ];
    for comp in comps do
        Q, mor := Explode(comp);
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
            recs := ReconstructionFromComponent(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
            Append(~recss, recs);
        end for;
        Append(~recsss, recss);
    end for;
    return recsss;
end if;

end intrinsic;


intrinsic IsogenyInformation(facss::.) -> .
{Returns homology exponents of isogeny induced by splitting, and tests if it is compatible with the various polarizations. For now, set AllIdems to true, AllPPs to false, ProjToIdem to true, and ProjToPP to true.}

Rs := &cat[ [ fac[2][2] : fac in facs ] : facs in facss ];
T := VerticalJoin(Rs);
S := SmithForm(ChangeRing(T, Integers()));
D := Diagonal(S);
exps := [ d : d in D | not d eq 1 ];

gYs := [ #Rows(R) div 2 : R in Rs ];
gX := &+gYs;
EYs := [ StandardSymplecticMatrix(gY) : gY in gYs ];
EX := StandardSymplecticMatrix(gX);

EY := DiagonalJoin(EYs);
test := IsMultiple(Transpose(T)*EY*T, EX);
return exps, test;

end intrinsic;
