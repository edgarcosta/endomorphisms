# Hyperelliptic tests:
R.<x> = PolynomialRing(QQ)
#F.<r> = NumberField(x^2 - 2)
#R.<x> = PolynomialRing(F)

# Curve input: specify g and h in its equation y^2 + h y = g.

Curves = [];

# Obvious CM
f = x^5-1
h = R(0)
C = HyperellipticCurve(f,h)
Curves.append(C)

# Product CM x Non CM
f = x^4-7
h = x^3+x
C = HyperellipticCurve(f,h)
Curves.append(C)

# Product CM x Non CM
f = 3*x^4+11*x^2+12
h = x^3
C = HyperellipticCurve(f,h)
Curves.append(C)

# Product of two non-CM elliptic curves
f = x^4+x^2
h = x^3+1
C = HyperellipticCurve(f,h)
Curves.append(C)

# Trivial
f = x^2+x
h = x^3+1
C = HyperellipticCurve(f,h)
Curves.append(C)

# Square of non-CM elliptic curve
f = x^5+x^4
h = x^3+x+1
C = HyperellipticCurve(f,h)
Curves.append(C)

# QM
f = -1
h = x^3+1
C = HyperellipticCurve(f,h)
Curves.append(C)


# Square of CM elliptic curve
f = x^5 -5*x^3-10*x^2-8*x-2
h = x^3
C = HyperellipticCurve(f,h)
Curves.append(C)


for c in Curves :
	End = EndomorphismData(c,100,molin_neurohr=true)
	print c;
	print "Numerically computed algebra:", End.geometric().algebra()
	print "Upper bound on the Z-rank of End(Jac(C)) over QQbar:", CurveRankBound(c), "\n\n"
