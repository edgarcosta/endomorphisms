AttachSpec("../endomorphisms/magma/spec");

CC := ComplexFieldExtra(500);
CCExtra := ComplexFieldExtra(1000);
D := -795;
Bs := BinaryQuadraticForms(D);
qs := ReducedForms(Bs);
qs := [ <1,1,199>, <3,3,67>, <5,5,41>, <15,15,17> ];
for q in qs do
    a, b, c := Explode(q);
    tau := (-b + Sqrt(CCExtra ! (b^2 - 4*a*c)))/(2*a);
    j := CC ! jInvariant([1, tau]);

    QQ := RationalsExtra();
    p := RelativeMinimalPolynomial(j, QQ);
    print p;

    RCC := PolynomialRing(CC);
    p := RCC ! p;
    print Minimum([ Abs(tup[1] - j) : tup in Roots(p, CC) ]);
end for;

exit;
