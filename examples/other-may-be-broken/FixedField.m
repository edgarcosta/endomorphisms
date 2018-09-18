/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 2);
prec := 100;
F := RationalsExtra(prec);
CC := F`CC;

f := x^48 - 8*x^46 + 38*x^44 - 258*x^42 + 1481*x^40 - 5354*x^38 + 13470*x^36 - 27892*x^34 + 52210*x^32 - 85828*x^30 + 120366*x^28 - 147390*x^26 + 158497*x^24 - 147390*x^22 + 120366*x^20 - 85828*x^18 + 52210*x^16 - 27892*x^14 + 13470*x^12 - 5354*x^10 + 1481*x^8 - 258*x^6 + 38*x^4 - 8*x^2 + 1;
R<x> := PolynomialRing(F);
K<r> := NumberFieldExtra(f);
print K`iota;
//0.7498517515247285760987411073793727808109426752799054404213560835339470492207108412109984363597047277 +
//0.6616058877725444755248069748760300689930335227136485754490972113117290556184778433681225959552561501*$.1

Gp, Gf, Gphi := AutomorphismGroupPari(K);
L := K;

Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
for H in Hs do
    print H;
    gens := [ Gphi(h) : h in Generators(H) ];
    K := FixedFieldExtra(L, gens);
    print K;
end for;

