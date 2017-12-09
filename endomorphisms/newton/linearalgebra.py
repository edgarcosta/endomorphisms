#  Copyright (C) 2015, 2016 Edgar Costa, Nicolas Mascot
#  See LICENSE file for license details.

from sage.all import ComplexField, Matrix, RealField, ZZ
from sage.all import floor, identity_matrix, random_vector, sqrt
from sage.rings.complex_field import is_ComplexField
def random_vector_in_span(M, n = None):
    m = M.nrows()
    if n is None:
        n = M.ncols()
    L = random_vector(M.base_ring() ,n);
    return [ sum( [ L[j]*M[i,j] for j in range(n) ] ) for i in range(m) ]

def float_sum(input):
    return float_sum2(sorted( input, key = lambda x: x.abs()))

def float_sum2(input):
    if len(input) <= 2:
        return sum(input)
    else:
        k  = floor(len(input)/2)
        return float_sum2(input[:k]) + float_sum2(input[k:])


def qrp(A, rank = None, extra_prec = 16):
    m = A.nrows()
    n = A.ncols()
    s = min(m,n)
    if rank is not None:
        #s = min(s, rank + 1);
        pass;
    base_prec = A.base_ring().precision();
    base_field = A.base_ring()
    if is_ComplexField(base_field):
        field = ComplexField(base_prec + extra_prec + m);
    else:
        field = RealField(base_prec + extra_prec + m);
#    field = base_field; 
    R = copy(A);
    A = A.change_ring(field);
    P = [ j for j in xrange(n) ];
    Q = identity_matrix(field,m)
#    print A.base_ring()
#    print Q.base_ring()
    
    normbound = field(2)**(-0.5*field.precision());

    colnorms2 = [ float_sum([ R[k,j].norm() for k in range(m)]) for j in range(n) ];
    for j in xrange(s):
        # find column with max norm
        maxnorm2 = colnorms2[j];
        p = j;
        for k in xrange(j, n):
            if colnorms2[k] > maxnorm2:
                p = k;
                maxnorm2 = colnorms2[k];     
        #it is worth it to recompute the maxnorm                
        xnorm2 = float_sum([ R[i,p].norm() for i in range(j, m)])                
        xnorm = sqrt(xnorm2);
        # swap j and p in A, P and colnorms
        if j != p:
            P[j], P[p] = P[p], P[j]
            colnorms2[j], colnorms2[p] = colnorms2[p], colnorms2[j]
            for i in xrange(0, m):
                R[i, j], R[i, p] = R[i, p], R[i, j];
        
        # compute householder vector
        u = [ R[k,j] for k in xrange(j, m) ];
        alpha = R[j,j];
        # alphanorm = alpha.abs();
        if alpha.real() < 0:
            u[0] -= xnorm
            z = -1
        else:
            u[0] += xnorm
            z = 1
            

        
        beta = xnorm*(xnorm + alpha * z); 
        if beta.abs() < normbound:
            beta = float_sum( [u[i].conjugate() * R[i+j,j] for i in range(m - j)]) 
        if beta == field(0):
            break;
        beta = 1/beta;
        
        # H = I - beta * u * u^*
        # Apply householder matrix
        #update R
        R[j,j] = -z * xnorm;
        for i in xrange(j+1,m):
            R[i,j] = 0
        for k in xrange(j+1, n):
            lambdR = beta * float_sum( [u[l - j].conjugate() * R[l, k] for l in xrange(j, m)] );
            for i in xrange(j, m):
                R[i, k] -= lambdR * u[i - j];
        
        #update Q
        for k in xrange(m):
            lambdQ = float_sum( [u[l - j].conjugate() * Q[l, k] for l in xrange(j, m)] );
            lambdQ *= beta
            for i in xrange(j,m):
                Q[i, k] -= lambdQ * u[i - j];
        
        for k in xrange(j + 1, n):
            square = R[j,k].norm();
            colnorms2[k] -= square;
            # sometimes is better to just recompute the norm
            if True: #colnorms2[k] < normbound or colnorms2[k] < square*normbound:
                colnorms2[k] = float_sum([ R[i,k].norm() for i in range(j + 1, m)])                
    # compute P matrix
    Pmatrix = Matrix(ZZ, n);
    for i in xrange(n):
        Pmatrix[ P[i], i] = 1;

    return Q.conjugate_transpose().change_ring(base_field), R.change_ring(base_field), Pmatrix;



def Kernel(A, rank):
    # used the QR factorization to extra dim elements of the kernel
    m = A.nrows();
    n = A.ncols();
    assert rank <= min(n,m)
    Q, R, P = qrp(A.conjugate_transpose(), rank);
    K = Matrix(A.base_ring(), n, n - rank)
    for i in xrange(K.nrows()):
        for j in xrange(K.ncols()):
            K[i, j] = Q[i, n - 1 - j]

    upper = R[rank - 1, rank - 1].abs();
    lower = 0;
    if n > rank and m > rank:
        lower = R[rank, rank].abs();
    return K, upper, lower;


def EqnMatrix(A, rank):
    m = A.nrows();
    n = A.ncols();
    assert rank <= min(n,m) 
    # Q = [Q1 Q2]
    # col(Q1) = col(A)
    # col(Q2) = col(A)^perp
    Q, R, P = qrp(A, rank);
    K = Matrix(A.base_ring(), m - rank, m);
    for i in xrange(m - rank):
        for j in xrange(m):
            K[ i, j] = Q[j, i + rank].conjugate();

    upper = R[rank - 1, rank - 1].abs();
    lower = 0;
    if n > rank and m > rank:
        lower = R[rank, rank].abs();
    return K,  upper, lower;



#qr with collumn pivoting for mpmath

#import sage.libs.mpmath.all as mpmath

import mpmath;

from copy import copy

from mpmath.libmp.backend import xrange


def qrp_mpmath(A, rank = None, extra_prec = 32, ctx = None):
    if ctx is None:
        ctx = A.ctx;
    # check values before continuing
    assert isinstance(A, ctx.matrix)
    m = A.rows
    n = A.cols
    assert n > 1
    assert m > 1
    assert extra_prec >= 0
    
    if rank is None:
        s = min(m, n);
    else:
        s = min(m, n, rank + 1);

    # check for complex data type
    cmplx = any(type(x) is ctx.mpc for x in A)

    normbound = mpmath.mpf(2)**(-0.75*mpmath.mp.prec);

    # temporarily increase the precision and initialize
    with ctx.extraprec(extra_prec):
        rzero = ctx.mpf('0.0')
        if cmplx:
            one = ctx.mpc('1.0', '0.0')
            def colnorm2(j, start = 0):
                if start > m:
                    return rzero
                else:
                    return ctx.re( ctx.fsum( R[i,j] * ctx.conj(R[i,j]) for i in xrange(start, m)) );
        else:
            one = ctx.mpf('1.0')
            def colnorm2(j, start = 0):
                if start > m:
                    return rzero
                else:
                    return  ctx.fsum( R[i,j]**2 for i in xrange(start, m) );
                    
        
        R = A.copy();
        Q = ctx.eye(m);
        P = [ j for j in xrange(0,n) ];  
            
        colnorms2 = [ colnorm2(j) for j in xrange(0, n) ];
        

        # main loop to factor A
        for j in xrange(0, s):
            
            
            maxnorm2 = colnorms2[j];
            p = j;
            for k in xrange(j, n):
                if colnorms2[k] > maxnorm2:
                    p = k;
                    maxnorm2 = colnorms2[k];     
                            
            xnorm2 = maxnorm2;
            xnorm = ctx.sqrt(xnorm2);
            
            
            
            # swap j and p in A, P and colnorms
            if j != p:
                P[j], P[p] = P[p], P[j]
                colnorms2[j], colnorms2[p] = colnorms2[p], colnorms2[j]
                for i in xrange(0, m):
                    R[i, j], R[i, p] = R[i, p], R[i, j];
                    
                    
            alpha = R[j,j]
            if cmplx:
                alphar = ctx.re(alpha)
                alphai = ctx.im(alpha)
                alphanorm2 = alphar**2 + alphai**2;
                alphanorm = ctx.sqrt(alphanorm2);
            else:
                alphanorm2 = alpha**2;
                if alpha > rzero:
                    alphanorm = alpha;
                else:
                    alphanorm = -alpha;
          
            # compute householder vector
            u = [ R[k,j] for k in xrange(j, m) ];
            
            
            z = one;
            if alphanorm != rzero:
                if cmplx:
                    z = alpha/alphanorm;
                else:
                    if alpha < rzero:
                        z = ctx.mpf("-1.0"); 
            u[0] += z * xnorm;

            
            beta = xnorm*(xnorm + alphanorm) ;
            if beta == 0:
                break;
            else:
                beta = 1/beta;
            
            # H = I - beta * u * u^*
            # Apply householder matrix
            #update R
            for k in xrange(j, n):
                if cmplx:
                    lambdR = beta * ctx.fsum( ctx.conj(u[l - j]) * R[l, k] for l in xrange(j, m) );
                else:
                    lambdR = beta * ctx.fsum( u[l - j] * R[l, k] for l in xrange(j, m) );

                for i in xrange(j, m):
                    R[i, k] -= lambdR * u[i - j];
                
            #update Q
            for k in xrange(m):
                if cmplx:
                    lambdQ = beta * ctx.fsum( ctx.conj(u[l - j]) * Q[l, k] for l in xrange(j, m) );
                else:
                    lambdQ =  beta * ctx.fsum( u[l - j] * Q[l, k] for l in xrange(j, m) )

                for i in xrange(j,m):
                    Q[i, k] -= lambdQ * u[i - j];
            
            for k in xrange(j + 1, n):
                if cmplx:
                    square =  ctx.re( R[j, k] * ctx.conj(R[j, k]) );
                else:
                    square = R[j,k]**2;
                colnorms2[k] -= square;
                if colnorms2[k] < normbound or colnorms2[k] < square*normbound:
                    colnorms2[k] = colnorm2(k, j + 1)
            
                        
        # compute P matrix
        Pmatrix = ctx.zeros(n);
        for i in xrange(n):
            Pmatrix[ P[i], i] = one;
        if cmplx:
            Qout = Q.transpose_conj();
        else:
            Qout = Q.transpose();
        return Qout, R, Pmatrix;




def qrp_wrapper(A, rank = None, extra_prec = 16, ctx = None):
    with mpmath.mp.workprec(A.base_ring().precision() + extra_prec):
#        print mpmath.mp.prec
        Am = mpmath.matrix(map(list, A.rows()))
        qm, rm, pm = qrp_mpmath(Am, rank, extra_prec, ctx)
        
        def to_sage(field,A):
            return Matrix(field, [[ field(A[i,j].real, A[i,j].imag) for j in range(A.cols)] for i in range(A.rows)])
        
        return to_sage(ComplexField(mpmath.mp.prec), qm), to_sage(A.base_ring(), rm), to_sage(A.base_ring(), pm)
#qrp = qrp_wrapper
