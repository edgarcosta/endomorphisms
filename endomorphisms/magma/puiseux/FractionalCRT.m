/***
 *  Chinese remainderings of fractions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


forward RandomSplitPrime;
forward RandomSplitPrimes;

forward FractionalCRTQQ;
forward FractionalCRTSplit;

forward ReduceConstantSplit;
forward ReducePointSplit;
forward ReduceMatrixSplit;
forward ReducePolynomialSplit;
forward ReduceRationalFunctionSplit;
forward ReduceBasisOfDifferentialsSplit;
forward ReduceAffinePatchSplit;
forward ReduceCurveSplit;


function RandomSplitPrime(f, B)
/*
 * Input:  A polynomial f over a general absolute field and a positive integer
 *         B
 * Output: A totally split prime of the rationals of size roughly 2^B and a
 *         root of f modulo that prime
 */

/* In the relative case things are more complicated, in the sense that we may
 * not actually get a totally reduced prime and the reduction below may fail.
 * The algorithm only sees the relative polynomial after all, which may not
 * give an absolute generator. Yet we need the current form since alternative
 * approaches like using IsTotallySplit seem to use much more time.
 *
 * One option, arguably the easiest, is to change base to an absolute field
 * before verifying. This is what we do for now. */

F := BaseRing(Parent(f));
while true do
    p := RandomPrime(B : Proof := false);
    FF := FiniteField(p);
    if Type(F) eq FldRat then
        test := true; rt := 1;
    else
        test, rt := HasRoot(MinimalPolynomial(F.1, Rationals()), FF);
    end if;
    if test then
        f_red := ReducePolynomialSplit(f, p, rt);
        if #Roots(f_red) eq Degree(f_red) then
            return [ p, Integers() ! rt ];
        end if;
    end if;
end while;

end function;


function RandomSplitPrimes(f, B, n)

ps_rts := [* *];
while #ps_rts lt n do
    p_rt_new := RandomSplitPrime(f, B);
    if not p_rt_new[1] in [ p_rt[1] : p_rt in ps_rts ] then
        Append(~ps_rts, p_rt_new);
    end if;
end while;
return ps_rts;

end function;


function FractionalCRTQQ(rs, prs)

rs := [ Integers() ! r : r in rs ];
M := Matrix([ [ 1, CRT(rs, prs), &*prs ] ]);
Lat := Lattice(Kernel(Transpose(M)));
v := ShortestVectors(Lat)[1];
return -v[1]/v[2];

end function;


function FractionalCRTSplit(rs, prs, OK, I, BOK, BI, K)
/*
 * rs is a set of remainders at the set of split primes prs,
 * I is an ideal in OK < K,
 * BOK and BI are eltseq of basis of OK and I
 */

if Type(OK) eq RngInt then
    return FractionalCRTQQ(rs, [ Norm(p) : p in prs ]);
end if;
n := CRT([ Integers() ! r : r in rs ], [ Norm(p) : p in prs ]);
M := Matrix(Integers(), [ [ b[i] : b in BOK ] cat [ KroneckerDelta(1, i)*n ] cat [ b[i] : b in BI ] : i in [1..#BOK] ]);
Lat := Lattice(Kernel(Transpose(M)));
v := ShortestVectors(Lat)[1];
if v[#BOK + 1] eq 0 then
    error "Division by zero";
end if;
return K ! ( (-1/v[#BOK + 1]) * &+[ v[i] * BOK[i] : i in [1..#BOK] ] );

end function;


function ReduceConstantSplit(x, p, rt)

FF := FiniteField(p); seq := Eltseq(x);
return &+[ (FF ! seq[i]) * rt^(i - 1) : i in [1..#seq] ];

end function;


function ReducePointSplit(P, p, rt);

return [ ReduceConstantSplit(c, p, rt) : c in Eltseq(P) ];

end function;


function ReduceMatrixSplit(M, p, rt);

return Matrix(FiniteField(p), [ [ ReduceConstantSplit(c, p, rt) : c in Eltseq(row) ] : row in Rows(M) ]);

end function;


function ReducePolynomialSplit(f, p, rt);

FF := FiniteField(p); R_red := PolynomialRing(FF, Rank(Parent(f)));
f_red := &+[ ReduceConstantSplit(MonomialCoefficient(f, mon), p, rt) * Monomial(R_red, Exponents(mon)) : mon in Monomials(f) ];
// Magma subtlety: we have to coerce to a univariate polynomial ring instead of
// a multivariate ring in one variable
if Rank(Parent(f)) eq 1 then
    return PolynomialRing(FF) ! f_red;
end if;
return f_red;

end function;


function ReduceRationalFunctionSplit(q, p, rt);

FF := FiniteField(p); R_red := PolynomialRing(FF, Rank(Parent(q)));
num_red := R_red ! ReducePolynomialSplit(Numerator(q), p, rt);
den_red := R_red ! ReducePolynomialSplit(Denominator(q), p, rt);
return num_red / den_red;

end function;


function ReduceBasisOfDifferentialsSplit(B, p, rt)

// Recall that we consider differentials as rational functions by fixing a
// uniformizer in a canonical way
FF := FiniteField(p); K_red := RationalFunctionField(FF, Rank(Parent(B[1])));
return [ K_red ! ReduceRationalFunctionSplit(b, p, rt) : b in B ];

end function;


function ReduceAffinePatchSplit(X, p, rt);

FF := FiniteField(p); R_red := PolynomialRing(FF, Rank(CoordinateRing(Ambient(X))));
DEs_red := [ R_red ! ReducePolynomialSplit(DE, p, rt) : DE in DefiningEquations(X) ];
return Curve(Scheme(AffineSpace(R_red), DEs_red));

end function;


function ReduceCurveSplit(X, p, rt)
/* NOTE: This only gives an affine patch (which is enough for our algorithms) */

U := ReduceAffinePatchSplit(X`U, p, rt); U`U := U;
U`is_hyperelliptic := X`is_hyperelliptic; U`is_planar := X`is_planar; U`is_smooth := X`is_smooth;
U`g := X`g; U`is_plane_quartic := X`is_plane_quartic;
U`P0 := U ! ReducePointSplit(X`P0, p, rt);
U`A := Ambient(U`U); U`RA := CoordinateRing(U`A); U`KA := FieldOfFractions(U`RA);
U`RU := CoordinateRing(U`U); U`KU := FunctionField(U`U); U`F := BaseRing(U`RU);
U`unif_index := X`unif_index;
U`DEs := DefiningEquations(U`U);
U`OurB := ReduceBasisOfDifferentialsSplit(X`OurB, p, rt);
U`NormB := ReduceBasisOfDifferentialsSplit(X`NormB, p, rt);
U`T := ReduceMatrixSplit(X`T, p, rt);
U`cantor_eqs := [* ReducePolynomialSplit(cantor_eq, p, rt) : cantor_eq in X`cantor_eqs *];
if assigned X`RRgens then
    nums := [ U`KU ! U`KA ! Numerator(ReduceRationalFunctionSplit(X`KA ! gen, p, rt)) : gen in X`RRgens ];
    dens := [ U`KU ! U`KA ! Denominator(ReduceRationalFunctionSplit(X`KA ! gen, p, rt)) : gen in X`RRgens ];
    U`RRgens := [ nums[i]/dens[i] : i in [1..#X`RRgens] ];
end if;
if assigned X`globgens then
    U`globgens := [ U`RA ! ReducePolynomialSplit(gen, p, rt) : gen in X`globgens ];
end if;
if assigned X`DEs_sub then
    U`DEs_sub := [ U`RA ! ReducePolynomialSplit(f, p, rt) : f in X`DEs_sub ];
end if;
U`initialized := true;
return U;

end function;
