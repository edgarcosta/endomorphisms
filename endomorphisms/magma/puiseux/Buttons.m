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

/*
// In principle, it would be better to do this without specifying P0, but 
// for a first go it's probably better to have the user provide this.
// 
// Doesn't work yet because in general we need to make the base change to the field
// over which endomorphisms are defined.
//
intrinsic CertifiedEndomorphismAlgebra(X::Crv, P0::Pt : Geometric := false, Al := "Divisor") -> .
{Returns the (certified) endomorphism algebra of X using the base point P0, 
by default over the base and over QQbar if Geometric is set to true. 
The output is the same as for HeuristicEndomorphismAlgebra with the last argument the certificates.
Al is either "Divisor" or "Cantor".}

  GeoEndoRep := GeometricEndomorphismRepresentation(X);
  if Geometric then
    EndoRep := GeoEndoRep;
  else
    F := BaseField(X);
    EndoRep := EndomorphismRepresentation(GeoEndoRep, F, hom<F -> F | >);
  end if;

  assert not IsWeierstrassPlace(Place(P0));

  certs := [* *];
  for MR in EndoRep do
    M := MR[1];
    if IsIdentity(MR[1]) then continue; end if; // no need to certify the identity
    if Al eq "Divisor" then
      bl, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M);
      Append(~certs, D);
    else
      bl, C := CantorFromMatrixAmbientSplit(X, P0, X, P0, M);
      Append(~certs, C);
    end if;
  end for;

  EndoAlg, EndoDesc := HeuristicEndomorphismAlgebra(X : Geometric := Geometric);
  
  return EndoAlg, EndoDesc, certs;
end intrinsic;
*/