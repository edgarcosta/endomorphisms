"""
 *  Bound functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

from sage.all import floor, PolynomialRing, QQ
import os
import inspect

def ComputeLPolys(C):
    return SmallJac(C)

#    f,h = C.hyperelliptic_polynomials()
#    d = discriminant(4*f+h^2)
#    LPolys = [0 for p in range(2,maxP)]
#    for p in range(2,maxP) :
#        if is_prime(p) and d%p != 0 :
#            Cp = C.base_extend(FiniteField(p))
#            LPolys[p] = Cp.frobenius_polynomial()
#       print p;
#    return LPolys

def SmallJac(C):
    f,h = C.hyperelliptic_polynomials()
    F = h^2+4*f
    g = floor((F.degree()-1) / 2)

    smalljacinvocation = str(F)
    smalljacinvocation = smalljacinvocation.replace(" ","")

    # TODO: Get everything pointing right and working
    filename = inspect.getframeinfo(inspect.currentframe())[0];
    __datadir__ = os.path.dirname(filename) + "/data/"
    if g >= 2:
        smalljacinvocation = __datadir__+"lpdata" + str(g) + " /tmp " + smalljacinvocation + " 256 >/dev/null";
    else :
        smalljacinvocation = __datadir__+"lpdata2 /tmp/bounds " + smalljacinvocation + " 256 >/dev/null";
    os.system(smalljacinvocation)
    with open("/tmp/bounds_lpdata.txt") as f:
        content = f.readlines()

    content = [x.strip() for x in content]
    R = PolynomialRing(QQ, 'x')

    content.remove(content[0])
    LPolys = [0 for p in range(2,257)]

    for l in content:
        if l.find("?") < 0 :
            LP = eval(l);
            p = LP[0];
            if g==3 :
                a3 = LP[1];
                a2 = LP[2];
                a1 = LP[3];
                LPolys[p] = x^6 + a3*x^5+a2*x^4+a1*x^3+p*a2*x^2+p^2*a3*x+p^3;

            if g==2 :
                a2 = LP[1];
                a1 = LP[2];
                LPolys[p] = x^4 + a2*x^3+a1*x^2+p*a2*x+p^2;

            if g==1 :
                a1 = LP[1];
                LPolys[p] = x^2 + a1*x+p;

    return LPolys;
