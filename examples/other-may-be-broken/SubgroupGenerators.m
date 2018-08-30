R<t> := PolynomialRing(Rationals());
f := t^2 - 2;
F<r> := NumberField(f);

R<t> := PolynomialRing(F);
L<s> := RelativeSplittingField(t^4 - (1 + 2*r));
print L;

Gp, Gf, Gphi := AutomorphismGroup(L);
print Gp.2;
K := FixedField(L, [ Gphi(Gp.2) ]);
print K;

print SubgroupGeneratorsUpToConjugacy(L, K, F);


F := Rationals();

R<t> := PolynomialRing(F);
L<s> := RelativeSplittingField(t^4 - 5);
print L;

Gp, Gf, Gphi := AutomorphismGroup(L);
print Gp.2;
K := FixedField(L, [ Gphi(Gp.2) ]);
print K;

print SubgroupGeneratorsUpToConjugacy(L, K, F);


exit;
