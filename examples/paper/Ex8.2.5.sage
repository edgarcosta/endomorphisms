"""
 *  An example from the paper
"""

from endomorphisms import EndomorphismData, mPlaneCurve

# Ambient ring:
F = QQ
P2.<x,y,z> = ProjectiveSpace(F, 2)

# Curve:
f = x^4 + 8*x^3*z + 2*x^2*y*z + 25*x^2*z^2 - x*y^3 + 2*x*y^2*z + 8*x*y*z^2 \
    + 36*x*z^3 + y^4 - 2*y^3*z + 5*y^2*z^2 + 9*y*z^3 + 20*z^4
X = mPlaneCurve(f)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)

print "Geometric representation:"
overK = Endo.over_field("geometric")
endodict = overK.full()
rep = overK.representation()
print rep

#print "Verifying endomorphisms:"
#test = Endo.verify()
#print test
