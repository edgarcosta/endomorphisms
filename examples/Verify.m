SetVerbose("EndoFind", 0);
SetVerbose("EndoCheck", 3);
SetVerbose("CurveRec", 0);

prec := 500;

// More examples over QQ
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f := -x^5;
h := x^3 + x + 1;
g := 4*f + h^2;
g *:= -1;
X := HyperellipticCurve(g);

f := x^6 + x^3 + 1;
f := x^6 + x^2 + 1;
f *:= -1;
X := HyperellipticCurve(f);

print "Verifying endomorphisms...";
test, fss := VerifyEndomorphismsLowerBound(X : Geometric := true);
print test;
