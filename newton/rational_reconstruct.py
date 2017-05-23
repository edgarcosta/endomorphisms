#  Copyright (C) 2016, Edgar Costa
#  See LICENSE file for license details

def rational_reconstruct_poly(p, m, n_deg = 0, d_deg = 0):
    """
    Construct a rational function n/d such that p*d is equivalent to n modulo m.

    returns (n, d)

    Examples:
        sage: P.<z> = QQ[]
        sage: p = -z**16 - z**15 - z**14 + z**13 + z**12 + z**11 - z**5 - z**4 - z**3 + z**2 + z + 1
        sage: m = z**21
        sage: n, d = rational_reconstruct_poly(p, m);
        sage: print n, d
        z**4 + 2*z**3 + 3*z**2 + 2*z + 1 z**10 + z**9 + z**8 + z**7 + z**6 + z**5 + z**4 +
        z**3 + z**2 + z + 1
        sage: print(p*d) % m == n
        True
    """
    # n and d are unique if m.degree() > (n.degree() + d.degree())
    if n_deg < 0 or d_deg < 0:
        raise ValueError("The degree bounds n_deg and d_deg should be positive.")

    if n_deg == 0:
        n_deg = (m.degree() - 1)//2
    if d_deg == 0:
        d_deg = (m.degree() - 1)//2

    P = p.parent()
    #XGCD until degree the degree of t1 surpasses the degree of n
    s0 = P(0);
    t0 = P(1);
    s1 = P(m);
    t1 = p % P(m);

    while n_deg < t1.degree():
        #not optimal, we should use divrem
        q = s1 // t1;
        r1 = s1 % t1;
        r0 = s0 - q*t0;
        s0 = t0
        s1 = t1
        t0 = r0
        t1 = r1

    assert(t0 != 0)
    if d_deg < t0.degree():
        raise ValueError("could not complete rational reconstruction")
    # make the denominator monic
    c = t0.leading_coefficient()
    t0 = t0.monic()
    t1 = t1/c
    return t1, t0
