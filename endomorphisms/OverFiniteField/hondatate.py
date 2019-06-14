#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.

# exposes some of Honda--Tate theory necessary to produce tight upper bounds
# For more details see Section 7.2

from utils import power_charpoly, tensor_charpoly
from sage.all import Integers, Rationals, LCM

def is_ordinary(f, q):
    r"""
    INPUT:

    - ``f`` -- the characteristic polynomial of frobenius

    - ``q`` -- the cardinality of the base field

    OUTPUT: True iff newton Polygon of ``f`` and the Hodge polygon match
    """

    if f(0) != 1:
        f = f.reverse();
    assert f(0) == 1

    flist = f.list();
    g = f.degree()//2;
    for elt in flist[:g]:
        if elt % q == 0:
            return False
    else:
        return True


def endomorphism_frob(f):
    r"""
    INPUT:

    -   ``f`` -- a Frobenius polynomial of an Abelian variety

    OUTPUT:

    A tuple:
            - dim_Q ( End(A^{al} )

            - k, the degree of field ext where all endomorphism are defined

            - a sorted list that represents the geometric isogeny decomposition
              over an algebraic closure.
              Let
                    A^{al} = (A_1)^n_1 x ... x (A_k)^n_t
              and write
                    det(1 - T Frob^k | H^1(Ai^n_i)) = c_i (T)^m_i.
              Then we return:
                    [ (m_i, m_i * deg(c_i), c_i) for i in range(1, t) ].

    EXAMPLE:

        sage: from endomorphisms.OverFiniteField import endomorphism_frob
        sage: ZZT.<T> = ZZ[]
        sage: F3 = 1 - T^2 + 9*T^4;
        sage: F7 = 1 + 4*T^2 + 49*T^4;
        sage: F13 = 1 - 8*T^2 + 169*T^4;
        sage: endomorphism_frob(F3)
        (8, 2, [(2, 4, 9*T^2 - T + 1)])
        sage: endomorphism_frob(F7)
        (8, 2, [(2, 4, 49*T^2 + 4*T + 1)])
        sage: endomorphism_frob(F13)
        (8, 2, [(2, 4, 169*T^2 - 8*T + 1)])

    """
    if f(0) != 1:
        f = f.reverse();
    assert f(0) == 1
    g = f.degree()/2;
    flist = f.list()
    q = Integers()(flist[-1]).nth_root(g);

    T = f.parent().gen()


    fof = tensor_charpoly(f,f);
    g = fof.change_ring(Rationals())(T/q)

    dimtotal = 0;
    fieldext = 1;

    for factor, power in g.factor():
        iscyclo = factor.is_cyclotomic(certificate=True)
        if iscyclo > 0:
            dimtotal += power *factor.degree()
            fieldext = LCM(fieldext, iscyclo)

    fext = power_charpoly(f, fieldext)

    endo = sorted([ (power, power * factor.degree(), factor) for factor, power in fext.factor()])

    return dimtotal, fieldext, endo
