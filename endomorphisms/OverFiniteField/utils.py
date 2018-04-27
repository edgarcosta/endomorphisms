#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.



from sage.all import PolynomialRing, CyclotomicField
# Functions regarding characteristic polynomials of linear transformations

def tensor_charpoly(f, g):
    r"""
    INPUT:

    - ``f`` -- the characteristic polynomial of a linear transformation

    - ``g`` -- the characteristic polynomial of a linear transformation

    OUTPUT: the characteristic polynomial of the tensor product of the linear transformations

    EXAMPLES::

    sage: x = PolynomialRing(ZZ,"x").gen();
    sage: tensor_charpoly((x - 3) * (x + 2),  (x - 7) * (x + 5))
    (x - 21) * (x - 10) * (x + 14) * (x + 15)

    """

    R = PolynomialRing(g.parent(), "y");
    y = R.gen();
    #x = g.parent().gen()
    A = f(y)
    B = R(g.homogenize(y))
    return B.resultant(A)

def power_charpoly(f, k):
    r""""
    INPUT:

    - ``f`` -- the characteristic polynomial of a linear transformation

    OUTPUT: the characteristic polynomial of the kth power of the linear transformation

    """
    K = CyclotomicField(k);
    R = f.base_ring();
    RT = f.parent();
    f = f.change_ring(K)
    T = f.parent().gen();
    a = K.gen()
    g = 1;
    for i in range(k):
        g *= f(a ** i  * T)
    glist = g.list();
    newf = [None]*(1 + f.degree());
    for i, ci in enumerate(glist):
        if i % k != 0:
            assert ci == 0, "i = %s, ci = %s" % (i, ci)
        else:
            newf[i/k] = R(ci)
    return RT(newf);

