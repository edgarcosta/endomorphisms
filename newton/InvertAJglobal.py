#  Copyright (C) 2015-2017 Edgar Costa
#  See the file LICENSE for license details.

import InvertAJlocal
import PicardGroup
import Divisor
import sage.libs.mpmath.all as mpmath
from sage.all import ComplexField, I
from sage.all import log

class InvertAJglobal:
    def __init__(self, f, prec, verbose = False):
        self.verbose = verbose;
        self.prec = prec;
        self.f = f;
        self.C = ComplexField(prec);
        self.Pic = PicardGroup.PicardGroup(f, self.C, self.verbose);
        self.iajlocal = InvertAJlocal.InvertAJlocal(f, prec, verbose = verbose);
        self.genus = self.iajlocal;
        
    def toC(self, x):
        return self.C(mpmath.mpmath_to_sage(x, self.prec));

    def toPoint(self, px):
        return ( self.toC(px) , self.toC( self.iajlocal.sign * mpmath.mp.sqrt( self.iajlocal.f(px) )) )


    def invertAJ(self, beta, k = 10):
        power2 = 2**k;
        beta_scaled = beta/power2;
        
        # the [1] coordinate is the error
        # the x coordinates of the inverse of the abel jacobi
        if self.verbose:
            print "Inverting locally with newton's method"
        
        with mpmath.extraprec(k*2 + 20):
            inverse_scaled_mpmath = self.iajlocal.solve(beta_scaled)[0];

        # convert x coordinates to points over self.C
        inverse_scaled_points = [ self.toPoint(px) for px in inverse_scaled_mpmath ];
        
        #base points
        iajbp = [ (self.toC(px), self.toC(py)) for px, py in self.iajlocal.get_basepoints()];
        if self.verbose:
            print "beta/2^power = Q - B where"
            print "Q = %s" % inverse_scaled_points;
            print "B = %s" % ([iajbp[0], iajbp[1]])
        
        # inverse_scaled_divisor = Q - B
        if self.verbose:
            print "Initializing divisors"
        inverse_scaled_divisor = Divisor.Divisor(self.Pic, [inverse_scaled_points[0], inverse_scaled_points[1] ]);
        basepoints_divisor = Divisor.Divisor(self.Pic, [ iajbp[0], iajbp[1] ]);

        if self.verbose:
            print "Computing their difference" 
        inverse_scaled_divisor -= basepoints_divisor;
        
        if self.verbose:
            print "Rescaling the difference"
        inverse = power2 * inverse_scaled_divisor

        return inverse;

    def test(self, iterations = 20, power = 10, randnorm = 1):

        for _ in range(iterations):
            try: 
                Qx = [ px + randnorm*(mpmath.rand() + I*mpmath.rand()) for px, _  in self.iajlocal.get_basepoints() ];

                QinJ = self.iajlocal.to_J_sum(Qx);

                Qtrace = Qx[0] + Qx[1];
                Qnorm = Qx[0] * Qx[1];

                Qiaj = self.invertAJ( QinJ, power);
                _, Qiajtrace, Qiajnorm = Qiaj.x_coordinate()
                Qiajtrace = -Qiajtrace;

                print "trace diff = %.3f norm diff = %.3f" % (mpmath.mp.prec + log(mpmath.fabs(Qtrace - Qiajtrace)/mpmath.fabs(Qtrace),2.0), mpmath.mp.prec + log(mpmath.fabs(Qiajnorm - Qnorm)/mpmath.fabs(Qnorm),2.0))
                if iterations == 1:
                    return max(mpmath.mp.prec + log(mpmath.fabs(Qtrace - Qiajtrace)/mpmath.fabs(Qtrace),2.0), mpmath.mp.prec + log(mpmath.fabs(Qiajnorm - Qnorm)/mpmath.fabs(Qnorm),2.0));
            except NotImplementedError as e:
                print e;
                print "Skipping iteration"
                pass



                    




        
