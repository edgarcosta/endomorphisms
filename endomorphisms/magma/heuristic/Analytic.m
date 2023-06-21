/***
 *  Finding an approximate basis for the geometric endomorphism ring
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

import "Recognition.m": MinimalPolynomialLLL;
forward GeometricEndomorphismRepresentationGH;
/* TODO: Algebraization is not completely a superstep of the complex calculation */

geo_hom_rep_CC := NewStore();

intrinsic CacheClearGeometricHomomorphismRepresentationCC()
{Clear the internal cache for GeometricHomomorphismRepresentationCC}
    // We need to save and restore the id, otherwise horrific things might
    // happen
    StoreClear(geo_hom_rep_CC);
    StoreSet(geo_hom_rep_CC, "cache", AssociativeArray());
end intrinsic;

intrinsic CacheGeometricHomomorphismRepresentationCC() -> Assoc
{Clear the internal cache for GeometricHomomorphismRepresentationCC}
    // We need to save and restore the id, otherwise horrific things might
    // happen
    if not StoreIsDefined(geo_hom_rep_CC, "cache") then
        StoreSet(geo_hom_rep_CC, "cache", AssociativeArray());
    end if;
    return StoreGet(geo_hom_rep_CC, "cache");
end intrinsic;

intrinsic ComplexStructure(P::ModMatFldElt) -> AlgMatElt
{Returns the complex structure that corresponds to the period matrix P. It is
the RR-linear transformation on homology that corresponds to multiplication
with i.}

CC := BaseRing(P); RR := RealField(CC);
PSplit := VerticalSplitMatrix(P); iPSplit := VerticalSplitMatrix(CC.1 * P);
return NumericalRightSolve(PSplit, iPSplit);

end intrinsic;


function RationalHomomorphismEquations(JP, JQ)
// Given two complex structures JP and JQ, returns the equations on homology
// satisfied by a homomorphism between the two corresponding abelian varieties.

/* Basic invariants */
RR := BaseRing(JP);
gP := #Rows(JP) div 2; gQ := #Rows(JQ) div 2;
n := 4 * gP * gQ;
/* Building a formal matrix corresponding to all possible integral
 * transformations of the lattice */
R := PolynomialRing(RR, n);
vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gQ, 2 * gP, vars);
/* Condition that integral transformation preserve the complex structure */
Comm := Eltseq(M * ChangeRing(JP, R) - ChangeRing(JQ, R) * M);
/* Splitting previous linear equations by formal variable */
return Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] : var in vars ]);

end function;


intrinsic TangentRepresentation(R::., P::ModMatFldElt, Q::ModMatFldElt : s0 := []) -> .
{Given a homology representation R of a homomorphism and two period matrices P
and Q, returns an analytic representation A of that same homomorphism, so that
A P = Q R. The keyword argument s0 is for repeated use with the same period
matrix P.}

CC := BaseRing(P);
if #s0 eq 0 then P0, s0 := InvertibleSubmatrix(P); end if;
P0 := Submatrix(P, [ 1..#Rows(P) ], s0);
QR := Q * ChangeRing(R, CC);
QR0 := Submatrix(QR, [ 1..#Rows(QR) ], s0);
A := NumericalLeftSolve(P0, QR0);
test := Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]);
if test gt CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;
return A, s0;

end intrinsic;


intrinsic TangentRepresentation(R::., P::ModMatFldElt : s0 := []) -> .
{Given a homology representation R of an endomorphism and a period matrix P,
returns an analytic representation of that same endomorphism, so that A P = Q
R.}

return TangentRepresentation(R, P, P : s0 := s0);

end intrinsic;


intrinsic HomologyRepresentation(A::., P::ModMatFldElt, Q::ModMatFldElt) -> .
{Given a complex tangent representation A of a homomorphism and two period
matrices P and Q, returns a homology representation R of that same
homomorphism, so that A P = Q R.}

CC := BaseRing(P);
SplitAP := VerticalSplitMatrix(A * P);
SplitQ := VerticalSplitMatrix(Q);
RRR := NumericalRightSolve(SplitQ, SplitAP);
R := Matrix(Rationals(), [ [ Round(cRR) : cRR in Eltseq(row) ] : row in Rows(RRR) ]);
test := Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]);
if test gt CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;
return R;

end intrinsic;


intrinsic HomologyRepresentation(A::., P::ModMatFldElt) -> .
{Given a complex tangent representation A of an endomorphism and a period
matrix P, returns a homology representation R of that same endomorphism, so
that A P = P R.}

return HomologyRepresentation(A, P, P);

end intrinsic;


intrinsic GeometricHomomorphismRepresentationCC(P::ModMatFldElt, Q::ModMatFldElt : s0 := []) -> SeqEnum
{Given period matrices P and Q, this function determines a ZZ-basis of
homomorphisms between the corresponding abelian varieties. These are returned
as pairs of a complex tangent representation A and a homology representation R
for which A P = Q R.}


/* Basic invariants */
gP := #Rows(P); gQ := #Rows(Q);

bool, cache := StoreIsDefined(geo_hom_rep_CC, "cache");
if not bool then
    cache := AssociativeArray();
end if;
precP := Precision(BaseRing(P));
precQ := Precision(BaseRing(Q));
k1 := <gP, precP, gQ, precQ>;
require Ncols(P) eq 2*gP: "P sould be a g x 2g matrix";
require Ncols(Q) eq 2*gQ: "Q sould be a g x 2g matrix";

if not IsDefined(cache, k1) then
    cache[k1] := AssociativeArray();
end if;

k2 := <P, Q>;
if not IsDefined(cache[k1], k2) then
    JP := ComplexStructure(P); JQ := ComplexStructure(Q);

    /* Determination of approximate endomorphisms by LLL */
    vprint EndoFind : "";
    vprint EndoFind : "Finding geometric homomorphisms...";
    M := RationalHomomorphismEquations(JP, JQ);
    Ker, test := IntegralLeftKernel(M);
    vprint EndoFind : "done finding geometric homomorphisms.";
    if not test then
        return [ ];
    end if;

    /* Deciding which rows to keep */
    CC := BaseRing(P); RR := BaseRing(JP); gens := [ ];
    for r in Rows(Ker) do
        R := Matrix(Rationals(), 2*gQ, 2*gP, Eltseq(r));
        /* Culling the correct transformations from holomorphy condition */
        Comm := ChangeRing(R, RR) * JP - JQ * ChangeRing(R, RR);
        if &and([ (RR ! Abs(c)) lt RR`epscomp : c in Eltseq(Comm) ]) then
            A, s0 := TangentRepresentation(R, P, Q : s0 := s0);
            Append(~gens, [* A, R *]);
        end if;
    end for;
    cache[k1][k2] := gens;
    StoreSet(geo_hom_rep_CC, "cache", cache);
end if;
return cache[k1][k2];

end intrinsic;


intrinsic GeometricEndomorphismRepresentationCC(P::ModMatFldElt : s0 := []) -> .
{Given a period matrix P, this function determines a ZZ-basis of endomorphisms
of the corresponding abelian variety. These are returned as pairs of a complex
tangent representation A and a homology representation R for which A P = P R.}

return GeometricHomomorphismRepresentationCC(P, P : s0 := []);

end intrinsic;


intrinsic GeometricHomomorphismRepresentation(P::ModMatFldElt, Q::ModMatFldElt, F::Fld : UpperBound := 16, DegreeDivides := Infinity(), s0 := []) -> SeqEnum
{Given period matrices P and Q, as well as a field F, this function determines
a ZZ-basis of homomorphisms between the corresponding abelian varieties. These
are returned as triples of an algebraized tangent representation A over a
number field K, a homology representation R and a complex tangent
representation ACC. We have ACC P = Q R, and via the infinite place of K the
matrix A is mapped to ACC. The inclusion of F into K is the second return
value.}

CC := BaseRing(P);
/* Determine matrices over CC */
gensPart := GeometricHomomorphismRepresentationCC(P, Q : s0 := []);
gensPart := [ [* ChangeRing(gen[1], F`CC), gen[2] *] : gen in gensPart ];
/* Determine minimal polynomials needed */
seqPart := &cat[ Eltseq(gen[1]) : gen in gensPart ];

vprint EndoFind : "";
vprint EndoFind : "Finding number field defined by homomorphisms...";
K, seq, hFK := NumberFieldExtra(seqPart, F : UpperBound := UpperBound, DegreeDivides := Infinity());
vprint EndoFind : "done finding number field defined by homomorphisms:", K;

assert #seq eq #seqPart;
if #seq eq 0 then
    return [ ], hFK;
end if;

r := #Rows(gensPart[1][1]); c := #Rows(Transpose(gensPart[1][1]));
As := [ Matrix(K, r, c, seq[((k - 1)*r*c + 1)..(k*r*c)]) : k in [1..#gensPart] ];
gens := [ [* As[k], gensPart[k][2] *] : k in [1..#gensPart] ];

/* Final check for correctness, removed because checking roots is more stable */
//for i in [1..#gens] do
//    abs := Max([ Abs(c) : c in Eltseq(EmbedMatrixExtra(gens[i][1]) - gensPart[i][1]) ]);
//    assert abs lt CC`epscomp;
//end for;
return gens, hFK;

end intrinsic;


intrinsic GeometricEndomorphismRepresentation(P::ModMatFldElt, F::Fld : UpperBound := 16, DegreeDivides := Infinity(), UseSplittingField := true, s0 := []) -> SeqEnum
{Given period matrices P and a field F, this function determines a ZZ-basis of
the corresponding abelian variety. These are returned as triples of an
algebraized tangent representation A over a number field K, a homology
representation R and a complex tangent representation ACC. We have ACC P = P R,
and via the infinite place of K the matrix A is mapped to ACC. The inclusion of
F into K is the second return value.}

Q := P;
/* Determine matrices over CC */
gensPart := GeometricHomomorphismRepresentationCC(P, Q : s0 := s0);
gensPart := [ [* ChangeRing(gen[1], F`CC), gen[2] *] : gen in gensPart ];
/* Determine minimal polynomials needed */
seqPart := &cat[ Eltseq(gen[1]) : gen in gensPart ];

/* Use splitting field instead of number field since the resulting field is
 * normal */
vprint EndoFind : "";
vprint EndoFind : "Finding number field defined by endomorphisms...";
if UseSplittingField then
    K, seq, hFK := SplittingFieldExtra(seqPart, F : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
else
    K, seq, hFK := NumberFieldExtra(seqPart, F : UpperBound := UpperBound, DegreeDivides := DegreeDivides);
end if;
vprint EndoFind : "done finding number field defined by endomorphisms:";
vprint EndoFind : K;

assert #seq eq #seqPart;
if #seq eq 0 then
    return [ ], hFK;
end if;

vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Minimal polynomials encountered during algebraization of tangent representations:";
vprint EndoFind, 3 : [ MinimalPolynomial(c) : c in seq ];

r := #Rows(gensPart[1][1]); c := #Rows(Transpose(gensPart[1][1]));
As := [ Matrix(K, r, c, seq[((k - 1)*r*c + 1)..(k*r*c)]) : k in [1..#gensPart] ];
gens := [ [* As[k], gensPart[k][2] *] : k in [1..#gensPart] ];

/* Final check for correctness, removed because checking roots is more stable */
//for i in [1..#gens] do
//    abs := Max([ Abs(c) : c in Eltseq(EmbedMatrixExtra(gens[i][1]) - gensPart[i][1]) ]);
//    assert abs lt BaseRing(P)`epscomp;
//end for;
return gens, hFK;

end intrinsic;


intrinsic GeometricEndomorphismRepresentationCC(X::. : UpperBound := 16, DegreeDivides := Infinity()) -> SeqEnum
{Given a curve X over a field F, this function determines a ZZ-basis of the
corresponding abelian variety. These are returned as triples of an algebraized
tangent representation A over a number field K, a homology representation R and
a complex tangent representation ACC. We have ACC P = P R for the period matrix
P of X, and via the infinite place of K the matrix A is mapped to ACC. The
inclusion of F into K is the second return value.}

assert ISA(Type(X), Crv);
P := PeriodMatrix(X); g := #Rows(P);
return GeometricEndomorphismRepresentationCC(P : s0 := [1..g]);

end intrinsic;


intrinsic GeometricEndomorphismRepresentation(X::. : UpperBound := 16, DegreeDivides := Infinity(), CC := false) -> SeqEnum
{Given a curve X over a field F, this function determines a ZZ-basis of the
corresponding abelian variety. These are returned as triples of an algebraized
tangent representation A over a number field K, a homology representation R and
a complex tangent representation ACC. We have ACC P = P R for the period matrix
P of X, and via the infinite place of K the matrix A is mapped to ACC. The
inclusion of F into K is the second return value.}

assert ISA(Type(X), Crv);
P := PeriodMatrix(X); g := #Rows(P); F := BaseRing(X);
if not CC and assigned X`geo_endo_rep then
    return X`geo_endo_rep;
end if;
if CC and assigned X`geo_endo_rep_CC then
    return X`geo_endo_rep_CC;
end if;

if CC then
    X`geo_endo_rep_CC := GeometricEndomorphismRepresentationCC(P : s0 := [1..g]);
    return X`geo_endo_rep_CC;
end if;
X`geo_endo_rep := GeometricEndomorphismRepresentation(P, F : UpperBound := UpperBound, DegreeDivides := DegreeDivides, s0 := [1..g]);
return X`geo_endo_rep;

end intrinsic;


/* Deprecated */
function GeometricEndomorphismRepresentationGH(q, f : CC := false)
/* Compability to calculate with geometrically hyperelliptic curves */

S := Parent(q); F := BaseRing(S); K := F;
S3 := S; K3 := FieldOfFractions(S3); PP2 := ProjectiveSpace(S3);
S2 := PolynomialRing(F, 2); K2 := FieldOfFractions(S2); PP1 := ProjectiveSpace(S2);
R := PolynomialRing(F); h21 := hom< S2 -> R | [ R.1, 1 ] >;

/* Find point on conic and parametrize */
Q := Conic(PP2, q);
test, P := HasRationalPoint(Q);
if not test then
    p := Evaluate(q, [ R.1, 0, 1 ]); K := NumberFieldExtra(p);
    S := PolynomialRing(K, 3); h := hom< Parent(q) -> S | [ S.1, S.2, S.3 ] >;
    q := h(q); f := h(f);

    S3 := S; KS3 := FieldOfFractions(S3); PP2 := ProjectiveSpace(S3);
    S2 := PolynomialRing(K, 2); KS2 := FieldOfFractions(S2); PP1 := ProjectiveSpace(S2);
    R := PolynomialRing(K); h21 := hom< S2 -> R | [ R.1, 1 ] >;

    Q := Conic(PP2, q); test, P := HasRationalPoint(Q);
end if;
phi := Parametrization(Q, P);
DE := DefiningEquations(phi);
h := hom< Parent(DE[1]) -> S2 | [ S2.1, S2.2 ] >;
DE := [ h(c) : c in DE ];

/* Create hyperelliptic curve */
// TODO: This could give incompatibilities! Not yet though.
g := Evaluate(f, DE);
PP2W := ProjectiveSpace(K, [Degree(f),1,1]);
S3W := CoordinateRing(PP2W);
h23 := hom< S2 -> S3W | [ S3.2, S3.3 ] >;
Y := Curve(PP2W, [ S3W.1^2 - h23(g) ]);

/* Create original curve */
PP3 := ProjectiveSpace(K, [2,1,1,1]);
S4 := CoordinateRing(PP3);
h34 := hom< S3 -> S4 | [ S4.2, S4.3, S4.4 ] >;
eqs := [ h34(q), S4.1^2 - h34(f) ];
X := Curve(PP3, eqs);

/* Isomorphism between the and corresponding pullback */
m := map< Y -> X | [ S3W.1, h23(DE[1]), h23(DE[2]), h23(DE[3]) ] >;
BX := BasisOfHolomorphicDifferentials(X);
BXPB := [ Pullback(m, omega) : omega in BX ];
BY := BasisOfHolomorphicDifferentials(Y);

/* Read off matrix on differentials */
genus := Degree(f) - 1;
rows := [ ];
for omega in BXPB do
    num := Numerator(omega/BY[#BY]);
    Rnum := Parent(num);
    row := [ MonomialCoefficient(num, Rnum.2^i) : i in [0..(genus - 1)] ];
    Append(~rows, row);
end for;
T := Matrix(rows); TCC := EmbedMatrixExtra(T);

/*
At this point, the period matrix for X is that of Y (taken wrt x^i dx / y), multiplied with T. So PX = T PY. Therefore if A PY = PY R, then
    A T^(-1) T PY = PY R
    T A T^(-1) T PY = T PY R
    T A T^(-1) PX = PX R
and T A T^(-1) is the new tangent representation.
*/

/* Hyperelliptic curve and its period matrix */
Y := HyperellipticCurve(h21(g)); PY := PeriodMatrix(Y);
PX := ChangeRing(TCC, BaseRing(Parent(PY))) * PY; g := #Rows(PX);

// TODO: Keep an eye on this
if CC then
    return GeometricEndomorphismRepresentationCC(PX : s0 := [1..g]), PX;
end if;
return GeometricEndomorphismRepresentation(PX, F : s0 := [1..g]), PX;

/* Possibly more stable: Find endomorphisms and take corresponding subfield */
GeoEndoRep := GeometricEndomorphismRepresentation(Y);
GeoEndoRep := [ [* T*tup[1]*T^(-1), tup[2] *] : tup in GeoEndoRep ];
LK := BaseRing(GeoEndoRep[1][1]);
seq := &cat[ Eltseq(tup[1]) : tup in GeoEndoRep ];
L, _, seq := SubfieldExtra(LK, seq);
GeoEndoRep := [ [* Matrix(genus, genus, seq[((i - 1)*genus^2 + 1)..(i*genus^2)]), GeoEndoRep[i][2] *] : i in [1..#GeoEndoRep] ];
return GeoEndoRep, PX;

end function;


intrinsic HeuristicEndomorphismFieldAsSplittingField(P::ModMatFldElt, F::Fld : UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns polynomials whose splitting field defines the field of definition of the endomorphisms.}

CC := BaseRing(P);
gensPart := GeometricEndomorphismRepresentationCC(P : s0 := []);
gensPart := [ [* ChangeRing(gen[1], F`CC), gen[2] *] : gen in gensPart ];
/* Determine minimal polynomials needed */
seqPart := &cat[ Eltseq(gen[1]) : gen in gensPart ];
prod := 1; minpols := [ ];
for i in [1..#seqPart] do
    //aCC := seqPart[Random([1..#seqPart])];
    aCC := seqPart[i];
    minpolnew := MinimalPolynomialLLL(aCC, F : DegreeDivides := DegreeDivides, UpperBound := UpperBound);
    prodnew := prod*minpolnew;
    if i eq 1 then
        prod := prodnew; G := GaloisGroup(prod); Append(~minpols, minpolnew);
        if #G eq DegreeDivides then return G, minpols; end if;
    else
        Gnew := GaloisGroup(prodnew);
        if #Gnew gt #G then
            G := Gnew; prod := prodnew; Append(~minpols, minpolnew);
            if #G eq DegreeDivides then return G, minpols; end if;
        end if;
    end if;
end for;
return G, minpols;

end intrinsic;


intrinsic HeuristicEndomorphismFieldAsSplittingField(X::Crv : UpperBound := 16, DegreeDivides := Infinity()) -> .
{Returns polynomials whose splitting field defines the field of definition of the endomorphisms.}

X := CurveExtra(X);
P := PeriodMatrix(X); F := BaseRing(X);
return HeuristicEndomorphismFieldAsSplittingField(P, F);

end intrinsic;
