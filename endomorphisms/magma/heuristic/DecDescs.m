/***
 *  Factor descriptions from lattices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function FactorDescriptionHyperelliptic(X, F)

if Genus(X) eq 1 then
    desc := "ell";

    K := BaseRing(X); F := BaseRing(K);
    K_seq := FieldDescription(K, F);
    field := K_seq;

    f, h := HyperellipticPolynomials(X);
    f_seq := Eltseq(f); h_seq := Eltseq(h);
    f_seq_seq := [ ElementDescription(coeff, F) : coeff in f_seq ];
    h_seq_seq := [ ElementDescription(coeff, F) : coeff in h_seq ];
    coeffs := [ f_seq_seq, h_seq_seq ];

else
    desc := "hyp";
    field := [ ];
    coeffs := [ ];
end if;

return [* desc, field, coeffs *];

end function;


function FactorDescriptionPlane(X, F)

desc := "pln";

K := BaseRing(X); F := BaseRing(K);

K_seq := FieldDescription(K, F);
field := K_seq;

f := DefiningPolynomials(X)[1];
mons := Monomials(f);
coeffsexps := [ ];
for mon in mons do
    coeff := MonomialCoefficient(f, mon);
    coeff_seq := ElementDescription(coeff, F);
    exp := Exponents(mon);
    Append(~coeffsexps, [* coeff_seq, exp *]);
end for;

return [* desc, field, coeffsexps *];

end function;


intrinsic FactorDescription(X::Crv, F::Fld) -> List
{Returns a string description of the curve X over the field F.}

if Type(X) eq CrvHyp then
    return FactorDescriptionHyperelliptic(X, F);
elif Type(X) eq CrvPln then
    return FactorDescriptionPlane(X, F);
else
    error "No implementation for general curves yet";
end if;

end intrinsic;
