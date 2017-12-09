R<x> := PolynomialRing(Rationals());
f := x^48 - 8*x^46 + 38*x^44 - 258*x^42 + 1481*x^40 - 5354*x^38 + 13470*x^36 - 27892*x^34 + 52210*x^32 - 85828*x^30 + 120366*x^28 - 147390*x^26 + 158497*x^24 - 147390*x^22 + 120366*x^20 - 85828*x^18 + 52210*x^16 - 27892*x^14 + 13470*x^12 - 5354*x^10 + 1481*x^8 - 258*x^6 + 38*x^4 - 8*x^2 + 1;
L<r> := NumberField(f);
DefineOrExtendInfinitePlace(L);
print L;

g := x^12 - 2*x^10 + 14/11*x^8 - 1/11*x^6 - 13/121*x^4 - 1/121*x^2 - 1/5324;

prec := 1000;
CC := ComplexFieldExtra(prec);
//SetEpsLLL(CC, 10^(-700));
aCC := Roots(g, CC)[1][1];
print aCC;

a := AlgebraizeElementInRelativeField(aCC, L);
print a;
print MinimalPolynomial(a);

exit;
