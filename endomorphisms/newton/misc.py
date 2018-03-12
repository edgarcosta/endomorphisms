#  Copyright (C) 2016-2017, Edgar Costa, Jeroen Sijsling
#  See LICENSE file for license details.

from sage.all import ComplexField, Infinity, Matrix, PolynomialRing, QQ, RR
from sage.all import algdep, floor, gp, is_square, sage_eval, sqrt, vector

def toCCap(x, prec = 53, ithcomplex_embedding = 0):
    try:
        return x.complex_embedding(prec = prec, i = ithcomplex_embedding);
    except AttributeError:
        return ComplexField(prec)(x);

def toCCap_list(x, prec = 53, ithcomplex_embedding = 0):
    return [toCCap(elt, prec, ithcomplex_embedding) for elt in x]

def almost_equal(z, z_alg, ithcomplex_embedding = 0, prec = None):
    if prec is None:
        prec = z.prec();

    z_algcc = toCCap(x = z_alg, prec = prec, ithcomplex_embedding = ithcomplex_embedding)

    if z_alg == 0:
        zabs = 1
    else:
        zabs = z_algcc.abs()
    return (z  - z_algcc).abs() < toCCap(2)**(-prec * 0.8) * zabs


def not_equal(z, z_alg, ithcomplex_embedding = 0, prec = None):
    if prec is None:
        prec = z.prec();

    z_algcc = toCCap(x = z_alg, prec = prec, ithcomplex_embedding = ithcomplex_embedding)

    if z_alg == 0:
        zabs = 1
    else:
        zabs = z_algcc.abs()

    return (z  - z_algcc).abs() > toCCap(2)**(-prec * 0.1) * zabs

def convert_magma_matrix_to_sage(m, K):
    try:
        return Matrix(K, m.sage());
    except:
        pass
    from re import sub
    text = str(m)
    text = sub('\[\s+', '[', text)
    text = sub('\s+', ' ', text)
    text = '['+text.replace(' ',', ').replace('\n',', ')+']'
    return Matrix(K, sage_eval(text, locals={'r': K.gens()}))

def NF_series_embedding(series, numberfield_series, ithcomplex_embedding = 0, known_bits=None ):
    nf = numberfield_series.base_ring()
    p = series.prec()
    out = [ ];
    for i, elt in  enumerate(series.list()):
        alg_elt = NF_embedding(elt, nf, ithcomplex_embedding, known_bits)
        if alg_elt is None:
            print "Truncating series at %s" % numberfield_series.gen()**i
            p = i;
            break;
        out += [alg_elt]

    return numberfield_series(out).add_bigoh(p)

def NF_embedding(z, numberfield, ithcomplex_embedding = 0, known_bits = None ):
    # NF_embedding_algdep is slower but more accurate
    return NF_embedding_algdep(z, numberfield, ithcomplex_embedding, known_bits)

def NF_embedding_algdep(z, numberfield, ithcomplex_embedding = 0, known_bits=None ):
    """
    Input:
        z \in Complex Number
        numberfield a NumberField
        ithcomplex_embedding uses the i-th embedding of self in the complex numbers
        known_bits to pass to algdep
    Output:
        returns z in terms of generators for the numberfield respecting the embedding

    Examples:
        QQx.<x> = QQ[]
        L.<a> = NumberField(x**2 + 3)
        NF_embedding_algdep(ComplexField(100)(2*sqrt(-3) + 1000000), L)
    """
    polynomialnf = PolynomialRing(numberfield, "W");

    if known_bits is None:
        prec = z.prec()*0.8
    else:
        prec = known_bits

    algdepnf = polynomialnf( algdep(z, numberfield.absolute_degree(), known_bits = known_bits) )
    #print algdepnf, algdepnf.roots()
    #if ithcomplex_embedding is None and numberfield.degree() > 1:
    #    ithcomplex_embedding = len(numberfield.complex_embeddings()) - 1

    minvalue = +Infinity;
    minarg = None;
    if numberfield.absolute_degree() > 1:
        for r, _ in algdepnf.roots():
            diff = (z  - r.complex_embedding(prec = prec, i = ithcomplex_embedding)).abs()
            if diff < minvalue:
                minvalue = diff;
                minarg = r
    else:
        minarg = algdepnf.roots()[0][0]
    assert almost_equal(z, minarg,  ithcomplex_embedding = ithcomplex_embedding, prec = known_bits), "z =? %s" % minarg;
    return minarg

def NF_embedding_linear(z, numberfield, ithcomplex_embedding = 0, known_bits = None ):
    """
    Input:
        z \in Complex Number
        numberfield a NumberField
        ithcomplex_embedding uses the i-th embedding of self in the complex numbers
        known_bits to pass to algdep
    Output:
        returns z in terms of generators for the numberfield respecting the embedding

    Examples:
        QQx.<x> = QQ[]
        L.<a> = NumberField(x**2 + 3)
        NF_embedding_algdep(ComplexField(100)(2*sqrt(-3) + 1000000), L)
    """

    if known_bits is None:
        prec = floor(z.prec()*0.8)
    else:
        prec = known_bits
    w = numberfield.gen();
    wcc = w.complex_embedding(prec = prec, i = ithcomplex_embedding);
    n = numberfield.degree();

    V=[z] + [wcc**i for i in range(n)]
    r = list(gp.lindep(V))
    print r
    if r[0] == 0:
        return None
    z_alg =  numberfield([QQ(-r[i]/r[0]) for i in range(1,n+1)])
    assert almost_equal(z, z_alg,  ithcomplex_embedding = ithcomplex_embedding, prec = known_bits);
    return z_alg



def find_rational_point(f, N = 2**8):
    out_rational = [];
    out_quad = [];
    for z in range(0, N):
        fz = f(z)
        fmz = f(-z)

        if is_square(fz) and fz != 0:
            out_rational += [(z, sqrt(fz))]
        elif fz != 0:
            out_quad += [(z, sqrt(fz))]
        if z != 0:
            if is_square(fmz) and fmz != 0:
                out_rational += [(-z, sqrt(fmz))]
            elif fmz != 0:
                out_quad += [(-z, sqrt(fmz))]

    return out_rational + sorted( out_quad, key=lambda point: RR(abs(point[1])) );


def branch(g, P):
    x0, y0 = P
    val = [(y0 - sqrt(g(x0))).norm(), (y0 + sqrt(g(x0))).norm() ];
    if val[0] < val[1]:
        return 1
    else:
        return -1



def sage_str_numberfield(nf, qq_variable, nf_variable):
    if nf == QQ:
        return 'QQ'
    else:
        assert nf.base_field() == QQ;
        Qpoly = PolynomialRing(QQ, qq_variable);
        return "NumberField(%s, '%s')" % (nf.defining_polynomial().subs(Qpoly.gen()), nf_variable);





def pade_linalg(input_series, m):
    ring = input_series.base_ring()
    polyring = input_series.parent()._poly_ring()
    an = vector(ring, input_series.padded_list())
    minus_an = -an;
    N = len(an) - 1
    n = N - m
    if n < 0:
        raise ValueError("the precision of the series is not large enough")
    Akj = Matrix(ring, N + 1, n + 1);
    for i in range(min(N + 1, n + 1)):
        Akj[i, i] = 1;

    Bkj = Matrix(ring, N + 1, m);
    for row in range(1, m+1):
        Bkj[row,:row] = minus_an[:row][::-1]
    for row in range(m+1, N+1):
        Bkj[row,:] = minus_an[row-m:row][::-1]
    C = Akj.augment(Bkj);
    pq = C.solve_right(an).list()
    p = pq[:n+1]
    q = [1] + pq[n+1:]
    return polyring(p), polyring(q)



def pade_sage(input_series, m, n):
    r"""
    Returns the Pade approximant of ``self`` of index `(m, n)`.

    The Pade approximant of index `(m, n)` of a formal power
    series `f` is the quotient `Q/P` of two polynomials `Q` and `P`
    such that `\deg(Q)\leq m`, `\deg(P)\leq n` and

    .. MATH::

        f(z) - Q(z)/P(z) = O(z**{m+n+1}).

    The formal power series `f` must be known up to order `n + m + 1`.

    See :wikipedia:`Pade\_approximant`

    INPUT:

    - ``m``, ``n`` -- integers, describing the degrees of the polynomials

    OUTPUT:

    a ratio of two polynomials

    .. WARNING::

        The current implementation uses a very slow algorithm and is not
        suitable for high orders.

    ALGORITHM:

    This method uses the formula as a quotient of two determinants.

    .. SEEALSO::

        * :mod:`sage.matrix.berlekamp_massey`,
        * :meth:`sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint.rational_reconstruct`

    EXAMPLES::

        sage: z = PowerSeriesRing(QQ, 'z').gen()
        sage: exp(z).pade(4, 0)
        1/24*z**4 + 1/6*z**3 + 1/2*z**2 + z + 1
        sage: exp(z).pade(1, 1)
        (-z - 2)/(z - 2)
        sage: exp(z).pade(3, 3)
        (-z**3 - 12*z**2 - 60*z - 120)/(z**3 - 12*z**2 + 60*z - 120)
        sage: log(1-z).pade(4, 4)
        (25/6*z**4 - 130/3*z**3 + 105*z**2 - 70*z)/(z**4 - 20*z**3 + 90*z**2
        - 140*z + 70)
        sage: sqrt(1+z).pade(3, 2)
        (1/6*z**3 + 3*z**2 + 8*z + 16/3)/(z**2 + 16/3*z + 16/3)
        sage: exp(2*z).pade(3, 3)
        (-z**3 - 6*z**2 - 15*z - 15)/(z**3 - 6*z**2 + 15*z - 15)

    TESTS:

    With real coefficients::

        sage: R.<z> = RR[[]]
        sage: f = exp(2*z)
        sage: f.pade(3, 3) # abs tol 1e-10
        (-1.0*z**3 - 6.0*z**2 - 15.0*z - 15.0)/(z**3 - 6.0*z**2 + 15.0*z - 15.0)

    When precision is too low::

        sage: f = z + O(z**6)
        sage: f.pade(4, 4)
        Traceback (most recent call last):
        ...
        ValueError: the precision of the series is not large enough
    """
    from sage.matrix.constructor import Matrix
    if input_series.precision_absolute() < n + m + 2:
        raise ValueError("the precision of the series is not large enough")
    polyring = input_series.parent()._poly_ring()
    z = polyring.gen()
    c = input_series.list()
    mat = Matrix(polyring, n + 1, n + 1)
    for i in range(1, n + 1):
        for j in range(n + 1):
            if m + i - j < len(c):
                mat[i, j] = c[m + i - j]
    for j in range(n + 1):
        mat[0, j] = z ** j
    resu_v = mat.determinant().truncate(n + 1)
    lead_v = resu_v.leading_coefficient()
    resu_v = resu_v / lead_v
    for j in range(n + 1):
        mat[0, j] = z ** j * (input_series.truncate(max(m - j + 1, 0)))
    resu_u = mat.determinant().truncate(m + 1)
    lead_u = resu_u.leading_coefficient()
    resu_u = resu_u / lead_u
    return lead_u / lead_v * resu_u / resu_v
