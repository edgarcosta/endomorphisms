#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.



from sage.all import GCD, Matrix, NumberField, PolynomialRing, PowerSeriesRing
from sage.all import gcd, sqrt
def sqrt_poly(g):
    if g.degree() == 0:
        return sqrt(g[0])

    d = GCD(g, g.derivative(g.parent().gen()))
    sg = g//d;
    return sg * sqrt_poly(g//sg**2)


def simplify_sqrt(c, d, rad):
    assert c.parent() == d.parent()
    assert rad.parent() == d.parent()
    # writes sqrt(c + d sqrt(rad))
    # as
    # a sqrt(d1) + b sqrt(d2)
    # where d1 * d2 = rad
    g2 = c**2 - d**2 * rad;
    g = sqrt_poly(g2)
    assert g**2 == g2, "%s != %s" % (g2, g**2)
    a = gcd( g + c, d)/2
    d1 = (g + c)/(2*a**2)
    b = d / (2 * a)
    d2 = rad/d1
    return a, b, d1, d2

def verify_algebraically_PS(g, P0, alpha, trace_and_norm, verbose = True):
    # input:
    # * P0 (only necessary to shift the series)
    # * [trace_numerator, trace_denominator, norm_numerator, norm_denominator]
    # output:
    # a boolean
    if verbose:
        print "verify_algebraically()"
    L = P0.base_ring()
    assert alpha.base_ring() is L
    L_poly = PolynomialRing(L, "xL");
    xL = L_poly.gen();
    # shifting the series makes our life easier
    trace_numerator, trace_denominator, norm_numerator, norm_denominator = [L_poly(coeff)(L_poly.gen() + P0[0]) for coeff in trace_and_norm];
    L_fpoly = L_poly.fraction_field()
    trace = L_fpoly(trace_numerator)/L_fpoly(trace_denominator)
    norm =  L_fpoly(norm_numerator)/L_fpoly(norm_denominator)
    
    Xpoly = L_poly([norm(0), -trace(0), 1])
    if verbose:
        print "xpoly = %s" % Xpoly

    if Xpoly.is_irreducible():
        M = Xpoly.root_field("c")
    else:
        # this avoids bifurcation later on in the code
        M = NumberField(xL, "c");
    if verbose:
        print M

    xi_degree = max([ elt.degree() for elt in [trace_denominator, norm_denominator] ])
    D = 2*xi_degree
    hard_bound = D + (4 + 2)
    soft_bound = hard_bound + 5
    M_ps = PowerSeriesRing(M, "T", default_prec = soft_bound);
    T = M_ps.gen();
    Tsub = T + P0[0];

    trace_M = M_ps(trace)
    norm_M = M_ps(norm)
    sqrtdisc = sqrt(trace_M**2 - 4 * norm_M)
    x1 = (trace_M - sqrtdisc) / 2
    x2 = (trace_M + sqrtdisc) / 2

    y1 = sqrt(g(x1))

    y2 = sqrt(g(x2))

    iy = 1/sqrt(g( Tsub ))

    dx1 = x1.derivative(T);
    dx2 = x2.derivative(T);

    dx1_y1 =  dx1/y1
    dx2_y2 = dx2/y2

    eq1 = Matrix([[-2*M_ps(alpha.row(0).list())(Tsub) * iy, dx1_y1, dx2_y2]])
    eq2 = Matrix([[-2*M_ps(alpha.row(1).list())(Tsub) * iy, x1 * dx1_y1, x2 * dx2_y2]])
    branches = Matrix([[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]] ).transpose()
    meq1 = eq1 * branches;
    meq2 = eq2 * branches;
    algzero = False
    for j in range(4):
        if meq1[0, j] == 0 and meq2[0, j] == 0:
            algzero = True;
            break;
    if verbose:
        print "Done, verify_algebraically()  = %s" % algzero
    return algzero

    
    
def verify_algebraically_GB(g, P0, alpha, trace_and_norm, verbose = True):
    # input:
    # * P0 (only necessary to shift the series)
    # * [trace_numerator, trace_denominator, norm_numerator, norm_denominator]
    # output:
    # a boolean
    if verbose:
        print "verify_algebraically()"
    L = P0.base_ring()
    assert alpha.base_ring() is L
    L_poly = PolynomialRing(L, "xL");
    xL = L_poly.gen();
    # shifting the series avoids makes our life easier
    trace_numerator, trace_denominator, norm_numerator, norm_denominator = [L_poly(coeff)(L_poly.gen() + P0[0]) for coeff in trace_and_norm];
    L_fpoly = L_poly.fraction_field()
    trace = L_fpoly(trace_numerator)/L_fpoly(trace_denominator)
    norm =  L_fpoly(norm_numerator)/L_fpoly(norm_denominator)
    L_fpoly_Z2 = PolynomialRing(L_fpoly, "z")
    z = L_fpoly_Z2.gen();
    gP0 = g(L_poly.gen() + P0[0])
    disc = (trace**2 - 4 * norm)


    IsqrtD = L_fpoly_Z2.ideal([z**2 - disc.numerator()])


    R0 = L_fpoly_Z2.quotient_ring(IsqrtD)

    iz = z/R0(z**2).lift();

    x1 = R0((trace - z / sqrt(disc.denominator()))/2).lift()
    x2 = R0((trace + z / sqrt(disc.denominator()))/2).lift()
    assert x1 + x2 == trace, "x1 + x2"
    assert R0(x1 * x2) == norm, "x1 * x2"

    dx1 = R0((trace.derivative(xL) - iz * sqrt(disc.denominator()) * disc.derivative(xL) /2)/2).lift()
    dx2 = R0((trace.derivative(xL) + iz * sqrt(disc.denominator()) * disc.derivative(xL) /2)/2).lift()
    assert R0(dx1 + dx2) ==  R0(trace.derivative(xL)), "dx1 + dx2"

    gx1 = R0(g(x1)).lift()
    gx2 = R0(g(x2)).lift()
    igx1 = R0( gx2 ).lift()/R0( gx2 * gx1 ).lift()
    assert R0(gx1 * igx1) == 1, "gx1 * igx1"
    igx2 = R0( gx1 ).lift()/R0( gx1 * gx2 ).lift()
    assert R0(gx2 * igx2) == 1, "gx2 * igx2"

    square = gx1.numerator()
    if verbose:
        print "Simplifying sqrt( g(x1).numerator() )"
    a1, a2, d1, d2 = simplify_sqrt(square.constant_coefficient().numerator(), square.monomial_coefficient(z).numerator(), disc.numerator())



    assert (d1.numerator()//gP0) in L or (d2.numerator()//gP0) in L, "d1 or d2"
    L_fpoly_Z = PolynomialRing(L_fpoly, 3, "z, y, w")
    z, y, w = L_fpoly_Z.gens();
    Isqrt = L_fpoly_Z.ideal([z - y * w,  y**2 - d1, w**2 - d2])
    R = L_fpoly_Z.quotient_ring(Isqrt)

    iz = z/(d1 * d2)
    assert R(z * iz) == 1
    iw = w/R(w**2).lift()
    assert R(w * iw) == 1
    iy = y/R(y**2).lift()
    assert R(iy * y) == 1


    sgx1 = R(a1 * y + a2 * w).lift()/R(sqrt(gx1.denominator())).lift()
    sgx2 = R(a1 * y - a2 * w).lift()/R(sqrt(gx1.denominator())).lift()
    assert R(sgx1**2) == gx1, "sgx1**2"
    assert R(sgx2**2) == gx2, "sgx2**2"

    isgx1 = R(sgx2).lift() / R(sgx1 * sgx2).lift()
    isgx2 = R(sgx1).lift() / R(sgx1 * sgx2).lift()

    assert R(sgx1 * isgx1) == 1, "sgx1 * isgx1"
    assert R(sgx2 * isgx2) == 1, "sgx2 * isgx2"
    
    if verbose:
        print "adjusting to d1//g(x) or d1/g(x) in L"

    if R(w**2/gP0).lift().degree() == 0:
        ct = iw * sqrt(R(w**2/gP0).lift().constant_coefficient())
    else:
        assert R(w**2/gP0).lift().degree() == 0
        ct = iy * sqrt(R(y**2/gP0).lift().constant_coefficient())

#    else:
#        assert R(y**2/gP0).lift() in L, "\n%s\n%s\n%s\n%s\n" % ( R(y**2/gP0).lift(),  R(y**2/gP0).lift() in L, R(w**2/gP0).lift(), R(w**2/gP0).lift() in L, )
#        ct = iy * sqrt(L(R(y**2/gP0).lift()))
    eq1 = Matrix([[-2*L_poly(alpha.row(0).list())(L_poly.gen() + P0[0]) * ct, R(dx1 * isgx1), R(dx2 * isgx2)]])
    eq2 = Matrix([[-2*L_poly(alpha.row(1).list())(L_poly.gen() + P0[0]) * ct, R(x1 * dx1 * isgx1), R(x2 * dx2 * isgx2)]])
    branches = Matrix(R, [[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, 1]] ).transpose()
    meq1 = eq1 * branches;
    meq2 = eq2 * branches;
    algzero = False
    for j in range(4):
        if meq1[0, j] == 0 and meq2[0, j] == 0:
            algzero = True;
            break;
    if verbose:
        print "Done, verify_algebraically()  = %s" % algzero
    return algzero


# TODO
# verify_algebraically in H^1 

def verify_algebraically(g, P0, alpha, trace_and_norm, verbose = True):
    return verify_algebraically_PS(g, P0, alpha, trace_and_norm, verbose)

