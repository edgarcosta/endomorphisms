function Lie(f, g);

return Derivative(f)*g - Derivative(g)*f;

end function;


function Det(G1, G2, G3);

G13 := Coefficient(G1, 2); G12 := Coefficient(G1, 1); G11 := Coefficient(G1, 0);
G23 := Coefficient(G2, 2); G22 := Coefficient(G2, 1); G21 := Coefficient(G2, 0);
G33 := Coefficient(G3, 2); G32 := Coefficient(G3, 1); G31 := Coefficient(G3, 0);
M := Matrix([
[G11, G12, G13],
[G21, G22, G23],
[G31, G32, G33]
]);
return Determinant(M);

end function;


function PI(G1, G2, G3);

return G1*G2*G3;

end function;


a := 2;
b := 3;
c := 5;
print 2^4*(c/a)*(b^2 - 4*a*c);

R<x> := PolynomialRing(Rationals());
F<r> := NumberField(x^3 - c);
F<r> := SplittingField((x^3 - c)*(x^6 - 2^4*(c/a)*(b^2 - 4*a*c)));
w := Roots(x^3 - c, F)[1][1];
R<x> := PolynomialRing(F);

f := (-x + w)/a;
g := (c - x^3)/a div f;
h := -b/a;

f := R ! f;
g := R ! g;
h := R ! h;

//R<x,y> := PolynomialRing(F, 2);
// y^3 = a x^4 + b x^2 + c
// a x^4 + b x^2 + (c - y^3) = 0
// y^4 + (b/a) y^2 + (c - x^3)/a = 0
//print 2*(y^4 - h*y^2 + f*g);

f2 := Coefficient(f, 2); f1 := Coefficient(f, 1); f0 := Coefficient(f, 0);
g2 := Coefficient(g, 2); g1 := Coefficient(g, 1); g0 := Coefficient(g, 0);
h2 := Coefficient(h, 2); h1 := Coefficient(h, 1); h0 := Coefficient(h, 0);

A := Matrix([ [ f2, f1, f0 ], [ h2, h1, h0 ], [ g2, g1, g0 ] ]);
B := A^(-1);

a1 := B[1,1]; b1 := B[1,2]; c1 := B[1,3];
a2 := B[2,1]; b2 := B[2,2]; c2 := B[2,3];
a3 := B[3,1]; b3 := B[3,2]; c3 := B[3,3];

a := a1 + 2*a2*x + a3*x^2;
b := b1 + 2*b2*x + b3*x^2;
c := c1 + 2*c2*x + c3*x^2;

f0 := 2*b*(b^2 - a*c);
X0 := HyperellipticCurve(f0);
I0 := G2Invariants(X0);
Fac := Factorization(f0);
print Fac;

F<r> := NumberField(Fac[2][1]);
print Subfields(F);

exit;
