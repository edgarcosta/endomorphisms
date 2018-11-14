/***
 *  Verifies that endomorphism ring is saturated
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic VerifySaturated(GeoEndoRep::SeqEnum, P::ModMatFldElt) -> BoolElt, AlgMatElt
{Returns a boolean that indicates whether the endomorphism ring in GeoEndoRep
is saturated in the corresponding algebra.}

Rs := [ gen[2] : gen in GeoEndoRep ];
// Creation of relevant algebras
g := #Rows(Rs[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring
A := Algebra(MatrixRing(Rationals(), 2*g)); GensA := [ A ! Eltseq(R) : R in Rs ];
// As a subalgebra
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
// As an associative algebra
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM); ind := Integers() ! ind;

ps := [ tup[1] : tup in Factorization(ind) ];
for p in ps do
    test, R := VerifySaturatedAtPrime(GeoEndoRep, P, p);
    if not test then
        return false, R;
    end if;
end for;
return true, "";

end intrinsic;


intrinsic VerifySaturatedAtPrime(GeoEndoRep::SeqEnum, P::ModMatFldElt, p::RngIntElt) -> BoolElt, AlgMatElt
{Returns a boolean that indicates whether the endomorphism ring in GeoEndoRep
is saturated in the corresponding algebra at p.}

/* Uses an extremely naive algorithm by excluding every intermediate element
 * directly */
Rs := [ gen[2] : gen in GeoEndoRep ];
CC := Parent(P[1,1]); RR := RealField(CC); JP := ComplexStructure(P);
I := [0..(p - 1)]; d := #Rs; CP := CartesianPower(I, d);

for tup in [ tup : tup in CP | not &and[ c eq 0 : c in tup ] ] do
    R := &+[ tup[i]*Rs[i] : i in [1..#Rs] ]/p;
    if &and[ c in Integers() : c in Eltseq(R) ] then
        D := JP*R - R*JP;
        if &and[ (RR ! Abs(d)) lt RR`epscomp : d in Eltseq(D) ] then
            return false, R;
        end if;
    end if;
end for;
return true, "";

end intrinsic;


intrinsic SaturateLattice(L::., M::. : ColumnsOrRows := "Columns") -> .
{Given a basis of a lattice L and a generating set of a lattice M in which L is
of finite index, returns a basis of M along with matrices that give expressions
of the provided generating sets in this basis. The flag ColumnsOrRows specifies
whether column of row vectors are interpreted as generating the lattice.}

/* In the end we have L = T B, M = U B in case of Rows and L = B T, M = B U in
 * case of Columns. */

if ColumnsOrRows eq "Columns" then
    B, T, U := SaturateLattice(Transpose(L), Transpose(M) : ColumnsOrRows := "Rows");
    return Transpose(B), Transpose(T), Transpose(U);
end if;

CC := BaseRing(L);
subL, s0 := InvertibleSubmatrix(L);
subM := Submatrix(M, [1..#Rows(M)], s0);
S := NumericalLeftSolve(subL, subM);
S, test := FractionalApproximationMatrix(S);
if not test then
    error "No suitable fractional approximation found";
end if;
/* At this point we have S L = M, where S has rational entries and an integral
* inverse */

/* Now we write S = R S0, where R is integral and where S0 has an integral
 * inverse */
S0 := ChangeRing(Matrix(Basis(Lattice(S))), Rationals());
S0i := S0^(-1); R := S * S0i;
/* The result is that the following B, T, U can be used */
B := ChangeRing(S0, CC) * L; T := S0i; U := R;

/* Final sanity check */
test1 := Minimum([ Abs(c) : c in Eltseq(L - ChangeRing(T, CC)*B) ]) lt CC`epscomp;
test2 := Minimum([ Abs(c) : c in Eltseq(M - ChangeRing(U, CC)*B) ]) lt CC`epscomp;
if not (test1 and test2) then
    error "Error in determining saturated lattice";
end if;
return B, T, U;

end intrinsic;
