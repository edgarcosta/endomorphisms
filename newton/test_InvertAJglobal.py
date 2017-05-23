#  Copyright (C) 2016-2017, Edgar Costa
#  See LICENSE file for license details

from sage.all import CC, Infinity, PolynomialRing, RealField, ZZ
from sage.all import ceil, copy, gcd, log, norm, set_random_seed, sqrt, vector
from AbelJacobi import AJ1, AJ2, PariIntNumInit, SetPariPrec
from Divisor import Divisor
from InvertAJglobal import InvertAJglobal

set_random_seed(1);
import random;
random.seed(1)
digits = 150;
prec = ceil(log(10)/log(2)*digits)
SetPariPrec(digits + 50);
PariIntNumInit(3);
power = 6
threshold = 2**(prec*.5)
ncurves = 8
#i = 8 needs more precision
for i in range(ncurves):
    print "i = %s of %s" % (i+1, ncurves)
    smoothQ = False;
    ZT = PolynomialRing(ZZ, "T");
    T = ZT.gen();

    while not smoothQ:
        g = ZT.random_element(degree = 6) # ZT([-4, 12, 0, -20, 0, 8, 1])# 
        dg = g.derivative()
        smoothQ = gcd(g, dg) == 1
    print "C: y^2 = %s" % g
    verbose = False 
    prec = ceil(log(10)/log(2)*digits)
    iaj = InvertAJglobal(g, prec, verbose)
    CCap = iaj.C;
    Pic = iaj.Pic
    
    y0 = 0
    while y0 == 0:
        x0 = 1#ZZ.random_element()
        y0 = sqrt(g(x0))
    print "\tPi = (%s, -/+%s)" % (x0, y0)
    P0, P1  =  vector(CCap, (x0, y0)), vector(CCap, (x0, -y0))

    y0 = 0
    while y0 == 0 or x0 == P0[0]:
        x0 = ZZ.random_element()
        y0 = sqrt(g(x0))
    print "\tQi = (%s, -/+%s)" % (x0, y0)
    Q0, Q1  =  vector(CCap, (x0, y0)), vector(CCap, (x0, -y0))
    
    points = [P0, P1, Q0, Q1]
    pairs = [[P0, Q0], [P0, Q1], [P1, Q0], [P1, Q1]]

    print "\tTesting P -> Divisor(P, {0,1} ) -> Divisor(P, {0,1}).compute_coordinates() = P - \inf_{(-1)^{0,1}}?"
    for i, P in enumerate(points):
        for sign in range(2):
            #point - \infty_{ (-1)^sign} ~
            # point + \inty_{ (-1)^{sign + 1} } - D0
            DP = Divisor(Pic, [P], sign);
            R0, R1 = DP.compute_coordinates()
            R0 = vector(R0)
            N = norm(P - R0)
            N = RealField(15)(N)
            print "\t\tP + \pm Infinity = R0 + \pm Infinity? : %s" % N
            if R1 != (-1)**(sign + 1) * Infinity or N > threshold:
                print "\t\t!! WARNING !!"
                print "\t\t\tsomething went wrong?"
                print "\t\t\tsign = %s R1 = %s" % (sign, R1)
                print "\t\t\tP = %s" % vector(CC, P)
                print "\t\t\tR0 = %s" % vector(CC, R0)
            
    print "\tTesting P, Q -> Divisor(P,Q) -> Divisor(P,Q).compute_coordinates() = P + Q - D0?"
    for i, pair in enumerate(pairs):
        p0, p1 = pair
        Dp0p1 = Divisor(Pic, [p0, p1]);
        R0, R1 = Dp0p1.coordinates();
        R0 = vector(R0)
        R1 = vector(R1)

        if norm(R0 - p0) > norm(R0 - p1):
            tmp = copy(R1)
            R1 = R0;
            R0 = tmp
        N = [norm(R0 - p0), norm(R1 - p1)]
        N = RealField(15)(max(N))

        print "\t\tP + Q = R0 + R1? :%s" % N
        if N > threshold:
            print "\t\t!! WARNING !!"
            print "\t\t\tsomething went wrong?"
            print "\t\t\tP = %s" % vector(CC, p0)
            print "\t\t\tQ = %s" % vector(CC, p1)
            print "\t\t\tR0 = %s" % vector(CC, R0)
            print "\t\t\tR1 = %s" % vector(CC, R1)

    print "\tTesting P, Q -> (Divisor(P,0)  + Divisor(Q,1)) -> (Divisor(P,0)  + Divisor(Q,1)).compute_coordinates() = P + Q - D0?"
    for i, pair in enumerate(pairs):
        p0, p1 = pair
        Dp0 = Divisor(Pic, [p0], 0)
        Dp1 = Divisor(Pic, [p1], 1);
        Dp0p1 = Dp0 + Dp1;
        R0, R1 = Dp0p1.coordinates();
        R0 = vector(R0)
        R1 = vector(R1)

        if norm(R0 - p0) > norm(R0 - p1):
            tmp = copy(R1)
            R1 = R0;
            R0 = tmp
        N = [norm(R0 - p0), norm(R1 - p1)]
        N = RealField(15)(max(N))

        print "\t\tP + Q = R0 + R1? :%s" % N
        if N > threshold:
            print "\t\t!! WARNING !!"
            print "\t\t\tsomething went wrong?"
            print "\t\t\tP = %s" % vector(CC, p0)
            print "\t\t\tQ = %s" % vector(CC, p1)
            print "\t\t\tR0 = %s" % vector(CC, R0)
            print "\t\t\tR1 = %s" % vector(CC, R1)

    print "\tTesting P -> (Divisor(P,0)  + Divisor(tau(P),1)) -> (Divisor(P,0)  + Divisor(tau(P),1)).compute_coordinates() = Q + tau(Q) - D0?"
    for i, pair in enumerate([[P0, P1], [Q0, Q1]] ):
        p0, p1 = pair
        Dp0 = Divisor(Pic, [p0], 0)
        Dp1 = Divisor(Pic, [p1], 1);
        Dp0p1 = Dp0 + Dp1;
        R0, R1 = Dp0p1.coordinates();
        if R1 != +Infinity and R1 != -Infinity:
            R0 = vector(R0)
            R1 = vector(R1)
            R1[1] = -R1[1]
        else:
            R1 = -R1

        if R0 == R1:
            print "\t\tP + tau(P) = Q + tau(Q)? True" 
        else:
            print "\t\tP + tau(P) = Q + tau(Q)? FALSE" 
            print "\t\t!! WARNING !!"
            print "\t\t\tsomething went wrong?"
            print "\t\t\tP = %s" % vector(CC, p0)
            print "\t\t\tQ = %s" % vector(CC, p1)
            print "\t\t\tR0 = %s" % (R0,)
            print "\t\t\ttau(R1) = %s" % (R1,)

    ajv = [ vector(CCap, AJ1(g, P)) for P in points ] 
    
    print "\tTesting P -> 2*aj(P) -> invertAJ(2*aj(P)) -> invertAJ(2*aj(P)).compute_coordinates() = P + P - D0?"
    for i, P in enumerate(points):
        double = iaj.invertAJ(2*ajv[i])
        R0, R1 = double.coordinates();
        R0 = vector(R0)
        R1 = vector(R1)
        N = [norm(R0 - R1), norm(R0 - P), norm(R1 - P)]
        N = RealField(15)(max(N))

        print "\t\t2*P = R0 + R1? :%s" % N
        if N > threshold:
            print "\t\t!! WARNING !!"
            print "\t\t\tsomething went wrong?"
            print "\t\t\tP = %s" % vector(CC, P)
            print "\t\t\tR0 = %s" % vector(CC, R0)
            print "\t\t\tR1 = %s" % vector(CC, R1)
    
    print "\tTesting P,Q -> aj(P,Q) -> invertAJ(aj(P) + aj(Q)) -> invertAJ(aj(P) + aj(Q)).compute_coordinates() = P + Q - D0?"
    for i, pair in enumerate(pairs):
        p0, p1 = pair
        aj = vector(CCap, AJ2(g, p0, p1))
        p0p1 = iaj.invertAJ(aj)
        R0, R1 = p0p1.coordinates();
        R0 = vector(R0)
        R1 = vector(R1)
        #print "R0 = %s" % vector(CC, R0)
        #print "R1 = %s" % vector(CC, R1)

        if norm(R0 - p0) > norm(R0 - p1):
            tmp = copy(R1)
            R1 = R0;
            R0 = tmp
        N = [norm(R0 - p0), norm(R1 - p1)]
        N = RealField(15)(max(N))

        print "\t\tP + Q = R0 + R1? :%s" % N
        if N > threshold:
            print "\t\t!! WARNING !!"
            print "\t\t\tsomething went wrong?"
            print "\t\t\tP = %s" % vector(CC, p0)
            print "\t\t\tQ = %s" % vector(CC, p1)
            print "\t\t\tR0 = %s" % vector(CC, R0)
            print "\t\t\tR1 = %s" % vector(CC, R1)
        



    





        
