#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

from sage.all import ComplexField, Infinity, Matrix, PolynomialRing, ZZ
from sage.all import copy, floor, log 

class Divisor:
    def __init__(self, *args, **kwargs):
        self.linearequations_cache = None;
        self.x_coordinate_cache = None;
        self.coordinates_cache = None;
        self.coordinates_input = False
        self.compute_coordinates_cache = None;
        self.sign = None;
        self.lower = 0;

        pic = args[0];

        if 'lower' in  kwargs:
            self.maxlower = kwargs['lower'];


        
        if 'basis_RR' in kwargs:
            self.init_basis_RR(pic, kwargs['basis_RR']);

        elif len(args) == 1:
            self.init_basis_RR(pic, pic.W0);

        elif len(args) == 2:
            oneargument = args[1];
            if type(oneargument) == type(Matrix(pic.R,[])):
                self.init_basis_RR(pic, oneargument)
            else:
                if len(oneargument) == 1:
                    self.init_point(pic, oneargument[0]);
                elif len(oneargument) == 2:
                    self.init_twopoints(args[0], oneargument[0], oneargument[1]);
                else:
                    raise NotImplemented;

        elif len(args) == 3:
            assert len(args[1]) == 1;
            oneargument = args[1];
            sign = args[2];
            assert len(oneargument) == 1;
            self.init_point(pic, oneargument[0], sign);            
        else:
            raise NotImplemented;

    def init_basis_RR(self, pic, basis_RR):
        self.Pic = pic;
        self.basis_RR = basis_RR;
        self.lower = max(self.lower, self.Pic.lower); 

    def init_point(self, pic, point, sign = 0):
        # Input: the picard group and a point in the Curve
        # Represents: the divisor class point - \infty_{ (-1)^sign }
        assert len(point) == 2;
        self.Pic = pic;
        self.coordinates_cache = [point];
        self.sign = sign;
        self.basis_RR, self.lower = self.Pic.PointMinusInfinity(point, sign);
        self.coordinates_input = True

    def init_twopoints(self, pic, P, Q):
        # Input: the picard group and two points in the Curve
        # Represents: the divisor class P + Q - \infty
        assert len(P) == len(Q);
        assert len(P) == 2;
        self.Pic = pic;
        self.coordinates_cache = [P, Q];
        self.basis_RR, self.lower = self.Pic.TwoPointsMinusInfinity(P, Q);
        self.coordinates_input = True

    


    def linearequations(self):
        # let A be the divisor represented by self
        # A is represented by H0(3D0 - D), more precisely, by its basis self.basis_RR in H(3D0)
        # where [D - D0] ~A
        # linearequations returns a basis for H0((g + 1)*\infty - E) in terms of H(D0)
        # where E is effective of degree g and 
        # E - (g/2) * \infty + D-(g+1)\infty ~ 0
        if self.linearequations_cache is None:
            self.linearequations_cache, lower  = self.Pic.linearequations(self.basis_RR.change_ring(self.Pic.Rextraprec));
            self.lower = max(self.lower, lower);
        return self.linearequations_cache;

    def x_coordinates(self):
        if self.x_coordinate_cache is None:
            c = self.coordinates()
            Rw =  PolynomialRing(self.Pic.Rdoubleextraprec, "T")
            G = Rw(1);
            for pair in c:
                if pair not in [+Infinity, -Infinity]:
                    G *= (Rw.gen() - pair[0])

            self.x_coordinate_cache = list(G)
        return self.x_coordinate_cache;

    def compute_coordinates(self):
        if self.compute_coordinates_cache is None:
            # recall that the zero divisor class A in Pic 
            # is represented by D such that [D - D0] ~ A
            # self.Pic.support(self.linearequations()) will return the support of E
            # where E is an effective divisor and
            # E - (g/2) * \infty + D - D0 ~ 0
            # E - (g/2) * \infty + A ~ 0
            # if A ~ [ P1 + ... + Pg - (g/2) \infty ]
            # thus E ~ g \infty  - (P1 + ... + Pg)
            # i.e E ~ tau(P1) + ... + tau(Pg) since P1 + tau(P1) ~ \infty 
            # tau = hyperelliptic involution
            supportE, lower = self.Pic.support(self.linearequations());
            self.lower = max(self.lower, lower);
            # now we need to apply tau, to recover P1, ..., Pg
            self.compute_coordinates_cache = [None]*self.Pic.g
            for i, P in enumerate(supportE):
                if P == +Infinity or P == -Infinity:
                    self.compute_coordinates_cache[i] = -P
                else:
                    assert len(P) == 2
                    self.compute_coordinates_cache[i] = (P[0], -P[1])
            
            
        return self.compute_coordinates_cache;
    
    def coordinates(self):
        if self.coordinates_cache is None:
            self.coordinates_cache = self.compute_coordinates();
        return self.coordinates_cache;

    

    def __repr__(self):
        out = "";
        out += "A divisor class of the hyperelliptic curve defined by %s" % self.Pic.f;
        if self.coordinates_input:
            if len(self.coordinates_cache) == 1:
                if self.sign == 0:
                    sign = "+"
                else:
                    sign = "-"
                out += " equivalent to [%s - \infty_%s]." % (self.coordinates_cache[0],sign);
            if len(self.coordinates_cache) == 2:
                 out += " equivalent to [%s + %s - (\infty_{+} + \infty_{-}]." % (self.coordinates_cache[0],self.coordinates_cache[1]);
        else:
            out += " equivalent to "
            for i in range(1, self.Pic.g):
                out += "P%d + " % i
            out += "P%d  - %d * (\infty_+ \infty_-) satisfying " % (self.Pic.g, self.Pic.g/2);
            out += "\nP_i = %s" % Matrix(ComplexField(20), self.compute_coordinates())
        return out

    
    def mDouble(self):
        tmp, lower = self.Pic.mDouble(self.basis_RR);
        return Divisor(self.Pic, basis_RR = tmp, lower = lower);

    def neg(self):
        tmp, lower =  self.Pic.Neg(self.basis_RR)
        return Divisor(self.Pic, basis_RR = tmp, lower = lower);

    def __neg__(self):
        return self.neg();

    def add(self, other):
        assert self.Pic == other.Pic;
        tmp, lower = self.Pic.Add(self.basis_RR, other.basis_RR)
        maxlower = max(lower, self.lower, other.lower);
        return Divisor(self.Pic, basis_RR = tmp, lower = maxlower );
    def __add__(self, other):
        return self.add(other);

    def sub(self, other):
        assert self.Pic == other.Pic;
        tmp, lower = self.Pic.Sub(self.basis_RR, other.basis_RR)
        maxlower = max(lower, self.lower, other.lower);
        return Divisor(self.Pic, basis_RR = tmp, lower = maxlower );

    def __sub__(self, other):
        return self.sub(other);



    
    def mul(self, a):
        assert a in ZZ;
        if a == 1:
            return copy(self)
        
        if a == 0:
            return Divisor(self.Pic);
        
        if a < 0:
            self.neg().mul( -a);
            a = -a;
            sign = 1;
        else:
            sign = 0;
        result = copy(self);
        negself = None;
        #negselflower = None;

        n = floor(log(a,2)) + 1;
        for i in range(n-2, -1, -1):
            result = result.mDouble();
            sign = (sign + 1)%2;
            if a & (1<<i):
                
                if sign == 1:
                    if negself == None:
                        negself = self.neg();

                    # result = -correct
                    tmp, lower =  self.Pic.AddFlip(result.basis_RR, negself.basis_RR);
                    lower = max(result.lower, lower, negself.lower);
                    result = Divisor(self.Pic, basis_RR = tmp, lower = lower);
                    #result = - (-correct - self) = self + correct
                    sign = 0;
                else:
                    #result = correct
                    tmp, lower = self.Pic.AddFlip(result.basis_RR, self.basis_RR);
                    lower = max(result.lower, lower, self.lower);
                    result = Divisor(self.Pic, basis_RR = tmp, lower = lower);

                    #result = -correct - self
                    sign = 1;
        
        if sign == 1:
            result = result.neg()
        
        return result
    
    
    
    def __mul__(a, b):
        if type(a) == type(2):
            return b.mul(a);
        elif type(b) == type(2):
            return a.mul(b);
        else:
            raise NotImplemented;

    def __rmul__(a, b):
        return a.mul(b);
