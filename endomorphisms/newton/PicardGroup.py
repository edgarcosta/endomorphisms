#  Copyright (C) 2015-2017  Edgar Costa, Nicolas Mascot
#  See LICENSE file for license details.

from linearalgebra import EqnMatrix, Kernel, qrp, random_vector_in_span
from sage.rings.complex_field import is_ComplexField as is_CF
from sage.all import ComplexField, Infinity, Matrix, PolynomialRing, RR, RealField, ZZ
from sage.all import ceil, copy, exp, identity_matrix, norm, sqrt, vector 

class PicardGroup:
    def __init__(self, f, ring, verbose = False):
        """
        input: f a polynomial over a number field
        ring:  Complex Field (with some prescribed precision) where we will work over


        We represent a divisor z \in Pic**0 by a divisor in D on X
        such that z = [D - D_0] where D_0 is a fixed divisor of degree d_0
        \deg D = d_0
        
        A divisor D on X is then represented by vector space H**0(\delta - D)
        for this representation to be faithful we need \delta - D to be base point free
        for that is sufficient \deg \delta >= \deg D + 2g, 
        as we will be dealing with divisors of degree d_0 and 2d_0 we pick \delta = 3 * D_0

        a vector space H**0(-) is represented by a bases of rational functions on X
        and these functions are represented by its evaluation at an effective divisor Z

        During our computations every function will lie in H**0(3D_0) or H**0(6D_0).
        Hence, for this representation to be faithful, it is necessary and sufficient that 
        H**0({3,6}D_0 - Z) = {0}, eg \deg D_0 = 6 d_0 + 1.
        If possible \supp D_0 \cap \supp Z = empty set would be ideal.

        In summary, a vector space W will be represented with \deg Z rows and \dim columns.

        """
        self.verbose = verbose;
        assert is_CF(ring)
        self.R = ring;
        self.prec = ring.precision()
        self.f = f;
        self.a0 = list(f)[0];
        self.an = list(f)[-1];
        self.g = int(ceil( f.degree() / 2 ) - 1); 
        self.d0 = 2 * self.g + 2;
        self.extraprec = self.prec + 3 * self.d0 + 16;
        self.Rextraprec = ComplexField(self.extraprec);
        self.Rdoubleextraprec =  ComplexField(2*self.extraprec)
        self.Pi = self.Rdoubleextraprec.pi();
        self.I = self.Rdoubleextraprec(0, 1);
        self.almostzero = self.R(2)**(- self.prec * 0.7);
        self.notzero = self.R(2)**(- self.prec * 0.1);
        
        self.lower = 0; 
        

        # D_0 = (g + 1)*( \infty_{-} + \infty_{+} ) = (g + 1) * \infty

        # Pick the effective divisor Z, ie, the evaluation points
        # Pick evaluation points /!\ Weierstrass pts !
        # in theory we only need \deg Z >= 6 d_0 + 1
        # in practice we aim for >= 6 d_0 + 2
        
        
        self.nZ = 6 * self.d0 + 2;
        if self.verbose:
            print "Using %s points on the unit circle, almost uniformly spaced, to represent each function" % ( self.nZ, )
        n = self.nZ;
        n = ceil(self.nZ/2);
        Z=[]
        
        while len(Z) < self.nZ:
            # X computed in self.Rdoubleextraprec
            X = [ exp( 2 * k * self.I * self.Pi / n ) for k in range(n) ]
            Z = [];
            for i, x in enumerate(X):
                #y = (-1)**i * sqrt( f(x) )
                #Z.append( (x, y) )
                #
                # x  in self.Rdoubleextraprec !
                y =  self.Rextraprec( sqrt( self.f(x) ) )
                x =  self.Rextraprec( x )
                # making sure the points don't collide
                if y.abs() < 1e-2:
                    Z.append( (x, y) )
                else:
                    Z.append( (x, y) )
                    Z.append( (x, -y) )
            n += 1
            

        self.Z = Z;
        self.nZ = len(Z); 
        #alternatively self.Z = Z[:self.nZ], but this might lose some of the good properties of having the points uniformly on the unit circle
            
        # V represents H**0 (3 D_0)
        V = Matrix(self.Rextraprec, self.nZ, 3 * self.d0 + 1 - self.g) # \dim = 6g + 6 + 1 - g = 5 g + 7

        # Vout represents H**0(D_0)
        Vout = Matrix(self.Rextraprec, self.nZ, self.d0 + 1 - self.g) # \dim = g + 3
        
        # H**0(D_0) = <1, x, x**2, ..., x**(g + 1), y>
        # H**0(3*D_0) = <1, x, x**2, ...,x**g, x**(g + 1), y, x**(g + 2) , y x, ..., x**(g + 1 + 2g + 2), y x**(2g + 2)>, dim = 5 g + 7 = 1 + g + 2*(2*g + 3)
        

        Vexps = [None]* (3 * self.d0 + 1 - self.g)
        i = 0;
        for n in xrange( 3 * (self.g + 1) + 1):
            Vexps[i] =  [ n, 0];
            i += 1;
            if n > self.g:
                Vexps[i] = [ n - self.g - 1, 1];
                i += 1;
        
        for i, (x, y) in enumerate(self.Z):
            for j in range( 3 * self.d0 + 1 - self.g):
                V[ i, j] = x ** Vexps[j][0] * y ** Vexps[j][1]

            for j in range( self.g + 2 ):
                Vout[i,j] = x**j

            Vout[ i, self.g + 2] = y

        self.V = V;
        self.Vout = Vout;
        self.Vexps = Vexps;
        
        # W0 represents H**0 (2 D_0)
        W0 = Matrix(self.Rextraprec, self.nZ, 2 * self.d0 + 1 - self.g) 
        for i in range( self.nZ ):
            for j in range( 2 * self.d0 + 1 - self.g):
                W0[i,j] = V[i,j]
        # KV = Ker(V***) = col(V)**perp 
        KV, upper, lower = EqnMatrix(V, 3 * self.d0 + 1 - self.g);
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        self.lower = max(self.lower, lower)

        self.KV = KV;
        self.W0 = W0;

    def threshold_check(self, upper, lower = None):
        # checking that upper != 0, and lower ~ 0 when compared with upper
        if self.verbose:
            print "upper  = %.10e lower = %.10e" % (upper, lower)
            if lower is not None:
                print (upper >= self.notzero) and (upper*self.almostzero >= lower)
        if lower is None:
            return upper >= self.notzero
        else:
            return (upper >= self.notzero) and (upper*self.almostzero >= lower);

    def AddFlip(self, A, B):

        maxlower = self.lower;

        field = self.R;
        output_field = self.R
        if A.base_ring() == self.Rextraprec or B.base_ring() == self.Rextraprec:
            field = self.Rextraprec;
        if A.base_ring() == self.Rextraprec and B.base_ring() == self.Rextraprec:
            output_field = self.Rextraprec;
        
        # Takes H0(3D0-A), H0(3D0-B), and returns H0(3D0-C), with C s.t. A+B+C~3D0
        # H0(6D0-A-B) = H0(3D0-A) * H0(3D0-B)

        if self.verbose:
            print "PicardGroup.AddFlip(self, A, B)"
        #### STEP 1 ####
        if self.verbose:
            print "Step 1";
        
        AB = Matrix(field, self.nZ, 4 * self.d0  + 1 - self.g);    
        KAB = Matrix(field, AB.nrows() - AB.ncols(), AB.nrows());
        
        ABncols = AB.ncols();
        ABnrows = AB.nrows();
        rank = AB.ncols();
 
        while True:
            # 5 picked at random
            l = 5;
            ABtmp =  Matrix(field, self.nZ, 4 * self.d0  + 1 - self.g + l);

            for j in range(ABtmp.ncols()):
                a = random_vector_in_span(A);
                b = random_vector_in_span(B);
                for i in range(ABnrows):
                    ABtmp[i, j] = a[i] * b[i];

            Q, R, P = qrp(ABtmp, ABncols );

            upper = R[rank - 1, rank - 1].abs();
            lower = 0;
            if ABtmp.ncols() > rank and ABnrows > rank:
                lower = R[rank, rank ].abs();
            
            if self.threshold_check(upper, lower):
                # the first AB.cols columns are a basis for H0(6D0-A-B) = H0(3D0-A) * H0(3D0-B)
                # v in AB iff KAB * v = 0
                for i in range(ABnrows):
                    for j in range(ABncols):
                        AB[i, j] = Q[i, j]
                    for j in range(ABnrows - rank):
                        KAB[j, i] = Q[i, j + rank].conjugate();
                maxlower = max(maxlower, lower);
                break;
            if self.verbose:
                print "Warning: In AddFlip, not full rank at step 1, retrying."
                print "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
                


        #### STEP 2 ####
        # H0(3D0-A-B) = {s in V : s*V c H0(6D0-A-B)}

        if self.verbose:
            print "Step 2"

        #### Prec Test ####
        if self.verbose:
            a = random_vector_in_span(A);
            b = random_vector_in_span(B);
            ab = Matrix(field, self.nZ, 1);
            for i in range(self.nZ):
                ab[i, 0] = a[i] * b[i];
            print field
            print "maxlower = %s ,lower = %s " % (RealField(15)(maxlower), RealField(15)(lower))
            print "Precision of the input of MakAddflip: %s"% (RealField(15)(max(map(lambda x: RR(x[0].abs()), list(KAB * ab)))),)

        KABnrows = KAB.nrows();
        KABncols = KAB.ncols();
        while True:
            # columns =  equations for H0(3D0-A-B) = {s in V : s*V c H0(6D0-A-B)} \subset  H0(6D0-A-B)
            l = 2;
            KABred = copy(KAB);
            KABred = KABred.stack( Matrix(field, l*KABnrows, KABncols ) );

            for m in range(l): # Attempt of IGS of V
                v = random_vector_in_span( self.V )
                for j in xrange(KABncols):
                    for i in range( KABnrows ):
                        KABred[i + m * KABnrows + KABnrows , j] = KAB[i, j] * v[j]

            ABred, upper, lower = Kernel(KABred, KABred.ncols() - (self.d0 + 1 - self.g));
            if self.threshold_check(upper, lower):
                maxlower = max(maxlower, lower);
                break;
            
            if self.verbose:
                print "Warning: In AddFlip, failed IGS at step 2, retrying."
                print "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))

        
        #### STEP 3 ####
        # H0(6D0-A-B-C) = f*H0(3D0), where f in H0(30-A-B)
        if self.verbose:
            print "Step 3"

        # fv represents f*H0(3D0)
        fV= copy(self.V)
        for j in xrange(3*self.d0 + 1 - self.g):
            for i in xrange(self.nZ):
                fV[i,j] *= ABred[i,0]

        KfV, upper, lower = EqnMatrix( fV, 3 * self.d0 + 1 - self.g ) # Equations for f*V
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = max(maxlower, lower);
 
        KfVnrows = KfV.nrows();
        KfVncols = KfV.ncols();

        ### STEP 4 ###
        # H0(3D0-C) = {s in V : s*H0(3D0-A-B) c H0(6D0-A-B-C)} = {s in V : s*H0(3D0-A-B) c f*V}
        if self.verbose:
            print "Step 4"
        while True:
            # Equations for H0(3D0-C)
            
            #expand KC
            l = 2;
            KC = self.KV.stack(Matrix(field, l * KfV.nrows(), self.KV.ncols()));
            KC_old_rows = self.KV.nrows();

            for m in range(l): # Attempt of IGS of H0(3D0-A-B)
                v = random_vector_in_span(ABred)
                for i in range(KfVnrows):
                    for j in range(KfVncols):
                        KC[i + m * KfVnrows + KC_old_rows, j] = v[j] * KfV[i, j];
            
            output, upper, lower = Kernel(KC, KC.ncols() - (2 * self.d0 + 1 - self.g) );
            if self.threshold_check(upper, lower):
                maxlower = max(maxlower, lower);
                break;

            if self.verbose:
                print "Warning: In MakAddFlip, failed IGS at step 4, retrying."
                print "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        
        return output.change_ring(output_field), maxlower;


    def Add(self, A, B):
        tmp, lower = self.AddFlip( A, B);
        tmp2, lower2 = self.AddFlip( tmp, self.W0);
        return tmp2, max(lower, lower2);
    
    def mDouble(self, A):
        return self.AddFlip( A, A);

    def Sub(self, A, B):
        tmp, lower =  self.AddFlip( A, self.W0)
        tmp2, lower2 =  self.AddFlip( tmp, B);
        return tmp2, max(lower, lower2);

    def Neg(self, A):
        return self.AddFlip( self.W0, A);

    def PointMinusInfinity(self, P, sign = 0):
        # Creates a divisor class of P- \inty_{(-1)**sign}
        maxlower = self.lower 
        assert self.f.degree() == 2*self.g + 2;
        sqrtan = (-1)**sign * self.Rextraprec( sqrt( self.Rdoubleextraprec(self.an)) );

        xP, yP = self.Rextraprec(P[0]),self.Rextraprec(P[1])
        assert (self.f(xP) - yP**2).abs() <= (yP**2).abs() * self.almostzero

        if xP != 0 and self.a0 != 0:
            x0 = self.Rextraprec(0);    
        else:
            x0 = xP
            while x0 == xP:
                x0 = self.Rextraprec(ZZ.random_element())
            
        y0 =  self.Rextraprec( sqrt( self.Rdoubleextraprec( self.f( x0 )))) # P0 = (x0, y0) 

        WP1 = Matrix(self.Rextraprec, self.nZ, 3 * self.g + 6) # H0(3D0 - P0 - g \infty) = H0(3D0 - P0 - g(\infty_{+} + \infty_{-}))
        EvP = Matrix(self.Rextraprec, 1, 3 * self.g + 6) # the functions of  H0(3D0-P0-g \infty) at P
        B =  Matrix(self.Rextraprec, self.nZ, 3 * self.g + 5) # H0(3D0 - \infty_{sign} - P0 - g\infty)

        # adds the functions x - x0, (x - x0)**1,  ..., (x - x0)**(2g+2) to it
        for j in xrange(2 * self.g + 2 ):
            for i, (x, y) in enumerate( self.Z ):
                WP1[ i, j] = B[ i, j] = (x - x0) ** (j + 1)
            EvP[ 0, j] = (xP - x0) ** (j + 1)
        
        # adds y - y0
        for i, (x, y) in enumerate( self.Z ):
            WP1[ i, 2 * self.g + 2 ] = B[ i, 2 * self.g + 2] = y - y0
        EvP[ 0, 2 * self.g + 2] = yP - y0

        # adds (x - x0) * y ... (x-x0)**(g+1)*y
        for j in range(1, self.g + 2):
            for i, (x,y) in enumerate( self.Z ):
                WP1[ i, 2 * self.g + 2 + j ] = B[ i, 2 * self.g + 2 + j] = (x-x0)**j * y 
            EvP[ 0, 2 * self.g + 2 + j] = (xP - x0)**j * yP


        # adds (x - x0)**(2g + 3) and  y *  (x - x0)**(g + 2)
        for i,(x,y) in enumerate(self.Z):
            WP1[ i, 3 * self.g + 4 ] =  (x - x0)**( 2 * self.g + 3)
            WP1[ i, 3 * self.g + 5 ] = y *  (x - x0)**(self.g + 2)
            B[ i, 3 * self.g + 4 ] =  (x - x0)**( self.g + 2) * ( y - sqrtan * x**(self.g + 1) )

        EvP[ 0, 3 * self.g + 4] = (xP - x0) ** ( 2 * self.g + 3)
        EvP[ 0, 3 * self.g + 5] = yP * (xP - x0)**(self.g + 2)
        
        # A = functions in  H0(3D0-P0-g \infty) that vanish at P
        # = H**0(3D0 - P0 - P -g \infty)
        K, upper, lower = Kernel(EvP, EvP.ncols() - (3 * self.g + 5))
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = max(maxlower, lower);

        A = WP1 * K
        # A - B = P0 + P + g \infty - (\infty_{+} + P0 + g\infty) = P - \inty_{+}
        res, lower = self.Sub(A,B)
        maxlower = max(maxlower, lower);
        return res.change_ring(self.R), maxlower;
    
    def TwoPointsMinusInfinity(self, P, Q): # Creates a divisor class of P + Q - g*\infty
        maxlower = self.lower
        xP, yP = self.Rextraprec(P[0]), self.Rextraprec(P[1])
        xQ, yQ = self.Rextraprec(Q[0]), self.Rextraprec(Q[1])

        assert (self.f(xP) - yP**2).abs() <= (yP**2).abs() * self.almostzero
        assert (self.f(xQ) - yQ**2).abs() <= (yQ**2).abs() * self.almostzero

        WPQ = Matrix(self.Rextraprec, self.nZ, 3 * self.g + 7) # H0(3D0  - g*\infty)
        Ev = Matrix(self.Rextraprec, 2, 3 * self.g + 7) # basis of  H0(3D0  - g*\infty) evaluated at P and Q
        Exps=[]
        for n in range( 2 * self.g + 4):
            Exps.append((n,0))
            if n > self.g:
                Exps.append(( n - self.g - 1, 1))

        for j in range( 3 * self.g + 7):
            u, v = Exps[j]
            for i in range( self.nZ ):
                x, y = self.Z[i]
                WPQ[ i, j ] = x**u * y**v

            Ev[0,j] = xP**u * yP**v
            Ev[1,j] = xQ**u * yQ**v
        
        K, upper, lower = Kernel(Ev, 2)
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = max(maxlower, lower);
        #Returns H0(3*D0 - P - Q - g infty)
        # note that  [P + Q + g infty - D0] ~ P + Q - infty
        return (WPQ * K).change_ring(self.R), maxlower


    def linearequations(self, W):
        maxlower = self.lower
        # From W = H0( 3*D0 - D ), return H0((g + 1)*\infty - E), where E is effective of degree g s.t. E - (g/2) * \infty + D-(g+1)\infty ~ 0. /!\ Assumes g even /!\
        # E0 = (3g/2+2)\infty = E + D + \infty
        # H0(3D0 - D - E0) = W( - E0) = {w in W : w*x**(3g/2+2) in W} has dim >= 1
        # phi \in W(-E0)  
        # <=>  div phi + 3D0 - D - E0  >= 0
        # <=>  div phi + 3D0 - D - E0 = E effective
        # <=> div phi = E + D - (3g/2 + 1)\infty
        assert self.g % 2 == 0, "the genus must be even for this implementation to work";

        KW, upper, lower = EqnMatrix(W.change_ring(self.Rextraprec), 2 * self.d0 + 1 - self.g);
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower));
        maxlower = max(maxlower, lower);

        KWnrows = KW.nrows();
        KWncols = KW.ncols();


        K = KW.stack(Matrix(self.Rextraprec, KWnrows, KWncols));

        for j, (x, y) in enumerate(self.Z):
            xpower =  x**( 3 * self.g/2 + 2);
            for i in range(KWnrows):
                K[i + KWnrows,j] = KW[i,j] * xpower;

        phi, upper, lower = Kernel(K, KWncols - 1);
        # rank >= 1
        if self.verbose:
            print "phi upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))    
        assert self.threshold_check(1,lower), "upper = %s lower = %s" % (RealField(35)(1), RealField(35)(lower))
        maxlower = max(maxlower, lower);

        # phi*H0( (5g/2+3) \infty) = H0( (4g + 4) \infty - D - E)
        Exps=[]

        for i in range( 5 * self.g / 2 + 4 ):
            Exps.append((i,0))
            if i > self.g:
                Exps.append(( i - self.g - 1, 1))
        
        # H0((4g+4)\infty-D-E)
        H4 = Matrix(self.Rextraprec, self.nZ, 4 * self.g + 7)

        for j, (u,v) in enumerate(Exps):
            for i, (x, y) in enumerate(self.Z):
                H4[i,j] = phi[i, 0] * x**u * y**v
        


        K4, upper, lower = EqnMatrix(H4, 4 * self.g + 7 )
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = max(maxlower, lower);


        # H0((g+1)\infty - E) = {s in V : s*W c H0((4g+4)OO-D-E)}
        # \dim >= 3 = self.d0 + 1 - 2 * self.g (equality holds for g <= 3)
        if self.g >= 4:
            print "Warning: Riemann--Roch correction term ignored when computing \dim H**0((g+1)\infty - E)!!!"
        K4nrows = K4.nrows();
        K4ncols = K4.ncols();
        while True:
            l = 2;
            Kres = copy(self.KV);
            Kres_old_rows = Kres.nrows();
            Kres = Kres.stack(Matrix(self.R, l * K4nrows, K4ncols));
            
            for m in range(l): # Attempt of IGS of W
                v=random_vector_in_span(W)
                for i in range(K4nrows):
                    for j in range(K4ncols):
                        Kres[i + m*K4nrows + Kres_old_rows, j] = v[j] * K4[i, j]
            
            res, upper, lower = Kernel(Kres, Kres.ncols() - (self.d0 + 1 - 2 * self.g));
            if self.threshold_check(upper, lower):
                maxlower = max(maxlower, lower);
                break;
            
            if self.verbose:
                print "Warning: In MakOutput, failed IGS at step 3, retrying."
                print "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        

        

        # Expressing H0((g+1)\infty-E) in terms of H0(D0)

        S = Matrix(self.R, self.nZ, (self.g + 3) + (self.d0 + 1 - 2 * self.g) )
        for i in range( self.nZ ):
            for j in range( self.g + 3):
                S[i,j] = self.Vout[i,j]

            for j in range( self.d0 + 1 - 2 * self.g):
                S[i ,j + self.g + 3] = res[i,  j ]

        K, upper, lower = Kernel(S, S.ncols() - ( self.d0 + 1 - 2 * self.g)) # S.cols - 3
        assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = max(maxlower, lower);
        
        # dropping the last  (self.d0 + 1 - 2 * self.g)  rows
        # K.rows = self.g + 3 
        out = K.matrix_from_rows([i for i in range(self.g + 3)])
        # a basis in H**0( D_0) for H**0((g+1)\infty - E)
        return out, maxlower
#
#        # Expressing H0((g+1)\infty-E) in terms of H0(3*D0)
#        # we could make have it in terms of H0(D0), however the extra rows are helpful for sanity checks
#
#        S = Matrix(self.R, self.nZ, (5*self.g + 7) + (self.d0 + 1 - 2 * self.g) )
#        for i in range( self.nZ ):
#            for j in range( 5*self.g + 7):
#                S[i,j] = self.V[i,j]
#
#            for j in range( self.d0 + 1 - 2 * self.g):
#                S[i ,j + 5*self.g + 7] = res[i,  j ]
#
#        K, upper, lower = Kernel(S, S.ncols() - ( self.d0 + 1 - 2 * self.g)) # S.cols - 3
#        assert self.threshold_check(upper, lower), "upper = %.10e lower = %.10e" % (upper, lower)
#
#        # dropping the last  (self.d0 + 1 - 2 * self.g)  rows
#        out = K.matrix_from_rows([i for i in range(5*self.g + 7)])
#        # K.rows = 5*self.g + 7
#        # a basis of  H**0((g+1)\infty - E) in terms of the basis elements of H**0( 3*D_0)
#        return out
    
    def support(self, basis):
        # Input:  a basis in H**0( D_0) for H**0((g+1)\infty - E), where E is an effective divisor of degree g/2
        # Output: returns the support of E as a list of g points 
 
        # Computing the intersection of 
        # <1, x, x**2,..., x**g> with <basis>
        B = basis.augment( identity_matrix(self.g + 1).stack(Matrix(basis.nrows() - (self.g + 1), (self.g + 1))))
        KB, upper, lower  = Kernel(B, basis.ncols() + self.g )
        if self.verbose:
            print "KB upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
        maxlower = self.lower
        
        result = [None for _ in range(self.g)]
        
        # if dim of the intersection of <1, x, x**2,..., x**g> with <basis> is 1
        if self.threshold_check(upper, lower):
            # we have the expected rank, all the roots are simple, unless E = 2*P
            # as in genus = 2 we can't have E = tau(P) + P
            maxlower = max(maxlower, lower);

            Rw = PolynomialRing(self.Rdoubleextraprec, "T")
            poly = Rw( KB.column(0)[-(self.g + 1):].list() )
            
            if self.g == 2:
                a, b, c = list(poly)
                disc = (b**2 - 4*a*c).abs()
                #a2 = (a**2).abs()
                dpoly = poly.derivative()
                if self.threshold_check((a**2).abs(), disc):
                    xcoordinates = dpoly.complex_roots() + dpoly.complex_roots()
                    maxlower = max(maxlower, disc);
                else:
                    xcoordinates = poly.complex_roots()
            else:
                xcoordinates = poly.complex_roots()

            for i, x in enumerate(xcoordinates):
                    y0 = self.Rextraprec(sqrt(self.f(x)));
                    x0 = self.Rextraprec(x)
                    vp = vector(self.R, [x0**u * y0**v for u, v in self.Vexps[:basis.nrows()] ]);
                    vn = vector(self.R, [x0**u * (-y0)**v for u, v in self.Vexps[:basis.nrows()] ]);
                    np = norm( vp * basis);
                    nn = norm( vn * basis);
                    #is it a simple root?
                    if dpoly(x).abs() < self.notzero:
                        assert self.g == 2, "one needs to be extra careful for genus > 2" 
                        
                    if self.threshold_check(np, nn) and not self.threshold_check(nn, np): #np > nn ~ 0
                        result[i] = (x0, -y0);
                        maxlower = max(maxlower, nn);
                    elif  not self.threshold_check(np, nn) and self.threshold_check(nn, np): #nn > np ~ 0
                        result[i] = (x0, y0);
                        maxlower = max(maxlower, np);

                    else:
                        print "sign for y is not clear: neg = %s, pos = %s,\ny = %s + %s * I = 0?" % tuple(vector(RealField(15), (nn, np, y0.real(), y0.imag())));
                        if nn > np:
                            result[i] = (x0, y0);
                            maxlower = max(maxlower, np);

                        else:
                            result[i] = (x0, -y0);
                            maxlower = max(maxlower, nn);

#                    else:
#                        assert self.g == 2, "for now only genus 2"
#                        # for g = 2,  we can't have E = P + tau(P), as (x - x0) \in H**0((g+1)\infty - E)
#                        # thus the rank of B should be at most  basis.ncols() + 1
#                        for j, xj in enumerate(xcoordinates[i+1:]):
#                            if (x0 - xj).abs() < self.notzero:
#                                # we might have a double root, and at the moment this should only work for g == 2
#                                break;
#                        #j is the other root 
#                        if self.threshold_check(np, nn) and not self.threshold_check(nn, np):
#                            result[j] = result[i] = (x0, -y0);
#                        elif not self.threshold_check(np, nn) and self.threshold_check(nn, np):
#                            result[j] = result[i] = (x0, y0);
#                        else:
#                            #this covers the case that y0 ~ 0
#                            if self.f(x0) > self.almostzero:
#                                print "sign for y is not clear: neg = %.2e, pos = %.2e,\nassuming y = %.2e + %.2e * I = 0" % (nn, np, y0.real(), y0.imag());
#                            result[j] = result[i] = (x0, 0);
#                    
#                                
            return result, maxlower;
        else:
            if self.verbose:
                print "there is at least a root at infinity, i.e. \inf_{+,-} \in supp E"
                print (RealField(35)(upper), RealField(35)(lower))
            # there is at least a root at infinity, i.e. \inf_{+,-} \in supp E
            # recall that \inf_{+} + \inf_{-} ~ P + tau(P), and in this case we have booth roots at infinity 
            assert self.g == 2, "for now only genus 2"


            Ebasis, upper, lower = EqnMatrix(basis,  basis.ncols() )
            assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
            maxlower = max(maxlower, lower);

            B = basis.augment( identity_matrix(self.g ).stack(Matrix(basis.nrows() - (self.g ), (self.g ))))
            KB, upper, lower  = Kernel(B, basis.ncols() + self.g - 1 )

            # 1 \in <basis>
            v0 = vector(self.R, [0] * basis.nrows())
            v0[0] = 1
            norm_one = norm(Ebasis * v0)
      
            if not self.threshold_check(upper, lower):
                # dim <1, x> cap <basis> = 2

                #1 \in <basis>
                assert self.threshold_check(1, norm_one),  "upper = %s lower = %s norm_one = %s " % (RealField(35)(upper), RealField(35)(lower), RealField(35)(norm_one));
                maxlower = max(maxlower, norm_one);
                
                #x \in <basis>

                vx = vector(self.R, [0]*basis.nrows())
                vx[1] = 1
                lower = norm(Ebasis * vx)
                assert self.threshold_check(1, lower), "lower = %s" % (RealField(35)(lower),);
                maxlower = max(maxlower, lower);
                
                vx2 = vector(self.R, [0]*basis.nrows())
                vx2[2] = 1
                lower = norm(Ebasis * vx2)
                if self.threshold_check(1, lower):
                    # E = \inf_{+} + \inf_{-} 
                    # <1, x, x^2>  \in H**0((g+1)\infty - E)
                    # x \in H^0 => 0 + 2 \infty - E >= 0
                    # x^2 \in H^0 => 4 * P_0 + \infty - E >= 0 for all P0

                    maxlower = max(maxlower, lower);
                    return [+Infinity, -Infinity], maxlower

                # <1, x>  \in H**0((g+1)\infty - E)
                # but not x^2
                # supp(E) \subsetneq { \inf_{+}, \inf_{-}  }
                if self.f.degree() % 2 == 0:
                    # we need to figure out the sign at infinity
                    Maug = Matrix(basis.nrows() - self.g - 1, self.g + 1).stack(identity_matrix(self.g + 1))
                    B = basis.augment( Maug )
                    KB, upper, lower  = Kernel(B, basis.ncols() + self.g )

                    a, b, c = KB.column(0)[-(2 + 1):].list()
                    assert self.threshold_check(upper, lower), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))

                    tmp = (-b/c) / sqrt(self.f.list()[-1]);
                    assert self.threshold_check( tmp.real(), tmp.imag()), "upper = %s lower = %s" % (RealField(35)(upper), RealField(35)(lower))
                    infsign = round(tmp.real());
                    assert infsign in [1, -1]


                    return [infsign*Infinity, infsign*Infinity], maxlower


                return [+Infinity, +Infinity], maxlower
            else:
                maxlower = max(maxlower, lower);
                #1 \notin <basis>
                assert self.threshold_check(norm_one),  "upper = %s lower = %s norm_one = %s " % (RealField(35)(upper), RealField(35)(lower), RealField(35)(norm_one));
                        
                        
            Rw = PolynomialRing(self.Rdoubleextraprec, "T")
            xcoordinates = Rw(KB.column(0)[-(self.g):].list()).complex_roots()
            assert len(xcoordinates) == 1;
            x0 = xcoordinates[0];
            y0 = self.Rextraprec(sqrt(self.f(x0)));

                        
            
            # if dim <1, x, x**2> cap <basis> >= 2
            # then dim <1, x> cap <basis> >= 1

            
            
            
            # we must have <1, x> \subset <basis>

            # double check that (x - x0) and (x - x0)**2 are in <basis>
            
            # (x - x0)
            v1 = vector(self.R, [0]*basis.nrows())
            v1[0] = -x0
            v1[1] = 1
            lower = norm(Ebasis * v1)
            assert self.threshold_check(1, lower), "lower = %s" % (RealField(35)(lower),);
            maxlower = max(maxlower, lower);

            # (x - x0)**2
            v2 = vector(self.R, [0]*basis.nrows())
            v2[0] = x0**2
            v2[1] = - 2*x0
            v2[2] = 1
            lower = norm(Ebasis * v2)
            assert self.threshold_check(1, lower), "lower = %s" % (RealField(35)(lower),);

            maxlower = max(maxlower, lower);

            
            # (x - x0)**3
            v3 = vector(self.R, [0]*basis.nrows())
            v3[0] = -x0**3
            v3[1] = 3*x0**2
            v3[2] = -3*x0
            v3[3] = 1
            nc = norm(Ebasis * v3)
            
            
            vp = vector(self.R, [x0**i * y0**j for i,j in self.Vexps[:basis.nrows()]]);
            vn = vector(self.R, [x0**i * (-y0)**j for i,j in self.Vexps[:basis.nrows()]]);

            np = norm( vp * basis);
            nn = norm( vn * basis);
            

            # either E = P + \inf_{+/-} 
            # or E = P + tau(P)   with P != \inf_{*}, but in this case E ~ \inf_+ + \inf_-        

            if self.threshold_check(np, nn) and not self.threshold_check(nn, np): #np > nn ~ 0
                result[0] = (x0, -y0);
                maxlower = max(maxlower, nn);

            elif not self.threshold_check(np, nn) and self.threshold_check(nn, np): #nn > np ~ 0
                result[0]  = (x0, y0);
                maxlower = max(maxlower, np);

            else:
                if self.threshold_check(1,nc) and self.threshold_check(1,np) and self.threshold_check(1,nn):
                    result[0] = (x0, y0)
                    result[1] = (x0, -y0)
                    maxlower = max(maxlower, nc, np, nn);
                    # this is equivalent to:  return result, maxlower
                    return [+Infinity, -Infinity], maxlower
                if self.f(x0) > self.almostzero :
                    print "sign for y is not clear: neg = %s, pos = %s,\nassuming y = %s + %s * I = 0" % tuple( vector( RealField(15), (nn, np, y0.real(), y0.imag())));
                result[0] = (x0, 0);
            
            # E = P + \inf_{+/-} 
            # P = result[0]
            x0, y0 = result[0]



            # Now figure out what infinity is in the supp of E
            # by checking if   y - s*x**(g+1) - (y0 - s*x0**(g+1)) \in <basis>
            # for s = sign * sqrt(an)
            sqrtan =  self.Rextraprec( sqrt( self.an) );
            # vp -> no pole at inf_{+} -> inf_{+} \in supp E
            vp = vector(self.R, [0]*basis.nrows())
            # vn -> no pole at inf_{-} -> inf_{-} \in supp E
            vn = vector(self.R, [0]*basis.nrows())
            vp[0] = -(y0 - sqrtan*x0**(self.g + 1))
            vn[0] = -(y0 + sqrtan*x0**(self.g + 1))
            vp[self.g + 1] = -sqrtan
            vn[self.g + 1] = sqrtan
            vp[self.g + 2] = 1
            vn[self.g + 2] = 1
            
            np = norm( Ebasis * vp );
            nn = norm( Ebasis * vn );

            


            if self.threshold_check(np * nc, nn) and not self.threshold_check(nn * nc, np) and not self.threshold_check(nn * np, nc): #np * nc > nn ~ 0
                maxlower = max(maxlower, nn);
                result[1] = -Infinity ;
            elif not self.threshold_check(np * nc, nn) and self.threshold_check(nn * nc, np) and not self.threshold_check(nn * np, nc): #nn * nc > np ~ 0
                maxlower = max(maxlower, np);
                result[1] = +Infinity
            elif not self.threshold_check(np * nc, nn) and not self.threshold_check(nn * nc, np) and self.threshold_check(nn * np, nc): # np*nn > nc
                # this should have been rulled out before, but now it looks more likely to have P + tau(P)
                # instead of E = P + inf
                maxlower = max(maxlower, nc);
                #this is equivalent to: result[1] = (x0, -y0)
                result = [+Infinity, -Infinity]
            else:
                print "inf_{?} \in supp E not clear: neg = %s, pos = %s, conj = %s" % tuple(vector(RealField(15),(nn, np, nc)));
                mn = min(nn, np, nc)
                maxlower = max(maxlower, mn);
                if np == mn:
                    result[1] = +Infinity;
                elif nn == mn:
                    result[1] = -Infinity;
                elif nc == mn:
                    #this is equivalent to: result[1] = (x0, -y0)
                    result = [+Infinity, -Infinity]
                else:
                    assert False, "something went wrong"
            return result, maxlower

            
            

            

            
            


       

