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
    F, h := BaseFieldExtra(BaseRing(GeoEndoRep[1][1]));
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
    F, h := BaseFieldExtra(BaseRing(GeoEndoRep[1][1]));
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
return Center(A) eq A;

end intrinsic;


intrinsic HeuristicIsGL2Generalized(X::Crv) -> .
{Returns whether or not X is of GL_2-type in the generalized sense of Booker--Sijsling--Sutherland--Voight--Yasaki.}

A := HeuristicEndomorphismAlgebra(X)[1];
if not Dimension(A) eq Genus(X) then
    return false;
end if;
if not Center(A) eq A then
    return false;
end if;
return #CentralIdempotents(A) eq 1;

end intrinsic;


intrinsic HeuristicJacobianFactors(X::Crv) -> .
{Returns factors of the Jacobian of X over the smallest possible fields, together with maps to these factors.}

P := PeriodMatrix(X);
GeoEndoRep := GeometricEndomorphismRepresentation(X);

comps := SplitComponents(P, GeoEndoRep);
Ys := [ ];
for comp in comps do
    Q, mor := Explode(comp);
    recs := ReconstructionsFromComponent(P, Q, mor);
    Ys cat:= recs;
end for;
return Ys;

end intrinsic;
