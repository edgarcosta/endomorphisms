#  Copyright (C) 2015-2017 Edgar Costa
#  See the file LICENSE for license details.
#
# TODO
# mpmath might be slow, only use it for the integrals...
# dynamically adjust precision of mpmath and sage
import sage.libs.mpmath.all as mpmath
from sage.all import CC, ComplexField, EllipticCurve, I, PolynomialRing, QQ, QQbar
from sage.all import cputime, floor, sqrt, walltime

class InvertAJlocal:
    def __init__(self, f, prec, verbose = False, method = 'gauss-legendre'):
        self.prec = prec;
        self.verbose = verbose;
        mpmath.mp.prec = prec;
        vars = f.variables()
        assert len(vars) == 1, "f must be a one variable polynomial"
        # for mpmath [cn, ... ,c0] represents cn x^n + ... + c0  
        self.fpoly = f.list();
        self.fpoly.reverse();
        self.f = lambda s: mpmath.polyval(self.fpoly, s);
        self.f_fp = lambda s: mpmath.fp.polyval(self.fpoly, s);

        self.genus = floor((f.degree() - 1)/2);
        
        self.basepoints = None;
        self.sign = None; 
        self.diff = None;
        self.sign = 1;
        self.method =  method;
        #if self.verbose:
        #    print "Initializing quadrature method";
        #    time A = mpmath.quad(lambda x: x, [0,1], method = self.method);
        #else:
        #    A = mpmath.quad(lambda x: x, [0,1], method = self.method);

    def random_good_point(self):
        with mpmath.extradps(100):
            while True:
                px = mpmath.rand() +  mpmath.rand()*I;
                py = self.sign * mpmath.sqrt(self.f(px));
                argpy = mpmath.fabs(mpmath.arg(py));
                if argpy < 3 and mpmath.fabs(py) > 1 and mpmath.almosteq(py**2, self.f(px)):
        		return (px, py) 


    def set_basepoints(self, basepoints = None):
        if not basepoints is None:
            assert len(basepoints) == self.genus
            px0, py0 = basepoints[0];
            if mpmath.sqrt(self.f(px0)) == py0:
                self.sign = 1;
            else:
                self.sign = -1

            for px, py in basepoints:
                assert mpmath.almosteq(self.f(px), py**2), "(%.3e, %.3e) is not a point on the curve" % (px, py, )
                assert mpmath.almosteq(mpmath.sqrt(self.f(px)), self.sign*py), "(%.3e, %.3e) doesn't has the same sign as (%.3e, %.3e) a point on the curve" % (px, py, px0, py0, )

            self.basepoints = basepoints;

        else:
            self.basepoints = [None] * self.genus;
            for i in range(self.genus): 
                self.basepoints[i] = self.random_good_point();

    def get_basepoints(self):
        if self.basepoints == None:
            self.set_basepoints();
        return self.basepoints;
         
        
    def jacobian(self, px_j):
        # d/dx \int x ^bx   
        sign = -self.sign;
        py_j = [ sign * mpmath.sqrt( self.f(px) ) for px in px_j ];
        return mpmath.matrix( [[ px_j[j] **i / py_j [j] for j in range(self.genus)] for i in range(self.genus)])
        
        
    def solve(self, beta, max_iterations = 100, tolerance = None, x0 = None):
        # solve \sum_j (\int_bj ^{P_j} w_i)_i = beta
        # assumes beta is small

        if tolerance == None:
            tolerance = mpmath.mpf(2)**(-mpmath.mp.prec + 3);
        #print "tolerance = %s" % tolerance 
        beta = mpmath.matrix([b for b in beta])
        while True:
            try:
                basepoints = self.get_basepoints();

                if x0 is None:
                    x = mpmath.matrix([ b[0]  for b in basepoints ]);
                    value = mpmath.matrix([0 for _ in range(self.genus)]);
                else:
                    assert len(x0) == self.genus
                    x =  mpmath.matrix(x0);

                    if self.verbose:
                        c, w = cputime(), walltime()
                        value = self.to_J_sum(x);
                        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
                    else:
                        value = self.to_J_sum(x);
                
                previous_norm = 1;
                divg_bound = 3;
                divergingQ = 0;
                for k in range(max_iterations):
                    J = self.jacobian(x);

                    # J y = f(x) - beta
                    y = mpmath.lu_solve(J, value - beta);
                    x = x - y;
                    
                    

                    norm = mpmath.norm(y)/mpmath.norm(x)
                    if self.verbose:
                        print "k = %d norm = %.3e" % (k, norm);
                    if norm < tolerance:
                        return (x,y);

                    if norm/previous_norm > 1.5:
                        divergingQ += 1;
                    else:
                        divergingQ = 0;

                    previous_norm = norm;

                    if divergingQ > divg_bound:
                        raise RuntimeWarning("method failed: the error diverged");
                    if self.verbose:
                        c, w = cputime(), walltime()
                        value = self.to_J_sum(x);
                        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
                    else:
                         value = self.to_J_sum(x);
                else:
                    raise RuntimeWarning("method failed: max number of iterations reached");

            except NotImplementedError as e:
                print e;
                print "Generating new basepoints"
                self.set_basepoints();
        
        
            
    def to_coordinates(self, px):
        return (px, self.sign * mpmath.sqrt(self.f(px)))


    def crossing_branchQ(self, px, bx):
        if self.verbose:
            print "crossing_branchQ"
        try:
            bp = bx - px;
            # first try to find a sign change
            for N in range(1, 10):
                last_eval = self.f_fp(px).imag;
                for i in range(1, 2**N + 1):
                    current_eval =  self.f_fp(px  + i * bp /(2**N)).imag;
                    if current_eval * last_eval < 0:
                        return True
                    else:
                        last_eval = current_eval

            root = mpmath.fp.findroot(lambda s: (self.f(px  + s * bp)).imag, (0, 1), solver = "secant");
            if (mpmath.fp.mpf(0) < root) and (root < mpmath.fp.mpf(1)):
                return True;
#                if self.f(bx  + root * (px - bx)).real < 2**(-mpmath.mp.prec/3):
#                    return True;
#                else:
#                   False 
            else:
                    return False;
        except ValueError:
            return False;

    def to_J_sum(self, x):
        if self.verbose:
            print "to_J_sum"
        assert len(x) == self.genus;
        basepoints = self.get_basepoints();
        value = [0] * self.genus;
        sign = self.sign;

        for j in range(self.genus): 
            px = x[j];
            bx = basepoints[j][0];
            if self.crossing_branchQ(px, bx):
                raise NotImplementedError("Path crosses a possible branchcut");

            path = lambda s: px  + s * (bx - px);

            for i in range(self.genus):
                diff = lambda s: sign * (bx - px) * ( s ** i ) / mpmath.sqrt( self.f( s ) )
                pathdiff = lambda s: diff( path(s) );
                value[i] += mpmath.quad(pathdiff, [0,1], method = self.method)
            
        return mpmath.matrix(value);
    

    def to_J(self, px, bx):
        # assuming px close to bx and py as the same sign as inf
        # I dont want to worry about flipping the sign of sqrt()
        sign = self.sign;
                     
        path = lambda s: px  + s * (bx - px);
                
        self.crossing_branchQ(px, bx);

        signdpath = sign * (bx - px);
        integrals_errors = [None] * self.genus;

        if self.crossing_branchQ(bx, px):
                raise NotImplementedError("Path crosses a possible branchcut");

        for i in range(self.genus):    
            pathdiff = lambda s: signdpath * ( path(s) ** i ) / mpmath.sqrt( self.f( path(s) ) );
            integrals_errors[i] = mpmath.quad(pathdiff, [0,1], error = True, method = self.method)

        coordinates = mpmath.matrix( [ x[0] for x in integrals_errors ] );
        error = max(map(lambda x: float(mpmath.nstr(x)),[x[1] for x in integrals_errors]));
        return (coordinates, error)
        
def test_elliptic_curve(a,b,alpha = 7, angles = 12):
    from sage.all import cos, sin, pi 
    E = EllipticCurve(QQ,[a,b])
    rationalpoints = E.point_search(10)
    EQbar = EllipticCurve(QQbar,[a,b])
    infx, infy, _ = rationalpoints[0]
    inf = EQbar(infx, infy)
    print inf
    R = PolynomialRing(QQ, "w")
    w = R.gen();
    f = w**3 + a*w + b        
    iaj = InvertAJlocal(f,256, method = 'gauss-legendre')
    iaj.set_basepoints( [(mpmath.mpf(infx),mpmath.mpf(infy)) ] )
    for eps in [cos(theta*2*pi/angles) + I*sin(theta*2*pi/angles)for theta in range(0,angles)]:
        px = infx + eps
        py = iaj.sign * sqrt(-E.defining_polynomial().subs(x = px, y = 0, z =1))
        P = EQbar(px, py)
        try:
            (v, errorv) = iaj.to_J(px, infx); 
            qx = (alpha*P - (alpha - 1)*inf).xy()[0]
            qxeps = CC(qx) + (mpmath.rand() + I*mpmath.rand())/1000
            (t1, errort1) = iaj.solve(alpha*v, 100, iaj.error, [qxeps]);
            qxguess = t1[0,0]
            print mpmath.fabs(qxguess - qx) < 2**(-mpmath.mp.prec/2)
        except RuntimeWarning as detail:
            print detail


