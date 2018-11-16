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
{What it says on the tin.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
return BaseRing(GeoEndoRep[1][1]);

end intrinsic;


intrinsic HeuristicEndomorphismAlgebra(X::Crv : Geometric := false) -> .
{What it says on the tin.}

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
{What it says on the tin.}

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
{What it says on the tin.}

GeoEndoRep := GeometricEndomorphismRepresentation(X);
return EndomorphismLattice(GeoEndoRep);

end intrinsic;


intrinsic HeuristicIsGL2Ribet(X::Crv) -> .
{What it says on the tin.}

A := HeuristicEndomorphismAlgebra(X)[1];
if not Dimension(A) eq Genus(X) then
    return false;
end if;
return Center(A) eq A;

end intrinsic;


intrinsic HeuristicIsGL2Generalized(X::Crv) -> .
{What it says on the tin.}

A := HeuristicEndomorphismAlgebra(X)[1];
if not Dimension(A) eq Genus(X) then
    return false;
end if;
if not Center(A) eq A then
    return false;
end if;
return #CentralIdempotents(A) eq 1;

end intrinsic;


/*
HeuristicJacobianFactors
*/
