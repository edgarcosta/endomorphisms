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
/* For a relative Galois extension defined as the splitting field of a
 * polynomial f */

F := BaseRing(Parent(f));
while true do
    p := RandomPrime(B : Proof := false);
    FF := FiniteField(p);
    // TODO: Only absolute number fields can be dealt with in this speedup.
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
    test := true;
    for p_rt in ps_rts do
        if p_rt[1] eq p_rt_new[1] then
            test := false;
            break;
        end if;
    end for;
    if test then
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
/* Need distinct primes downstairs and do not check for this */

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

U := ReduceAffinePatchSplit(X`U, p, rt);
U`is_hyperelliptic := X`is_hyperelliptic; U`is_planar := X`is_planar; U`is_smooth := X`is_smooth;
U`g := X`g; U`is_plane_quartic := X`is_plane_quartic;
U`U := U;
U`P0 := U ! ReducePointSplit(X`P0, p, rt);
U`patch_index := X`patch_index;
U`A := Ambient(U`U); U`R := CoordinateRing(U`A);
if X`is_hyperelliptic and (X`patch_index eq 3) then
    U`x := U`R.2; U`y := U`R.1;
else
    U`x := U`R.1; U`y := U`R.2;
end if;
U`F := BaseRing(U`R); U`K := FieldOfFractions(U`R);
U`DEs := DefiningEquations(U`U);
U`unif_index := X`unif_index; U`unif := (U`R).(U`unif_index);
U`OurB := ReduceBasisOfDifferentialsSplit(X`OurB, p, rt);
U`NormB := ReduceBasisOfDifferentialsSplit(X`NormB, p, rt);
U`T := ReduceMatrixSplit(X`T, p, rt);
if U`is_planar then
    U`cantor_equations := [* ReducePolynomialSplit(cantor_eq, p, rt) : cantor_eq in X`cantor_equations *];
end if;
U`initialized := true;
return U;

end function;
