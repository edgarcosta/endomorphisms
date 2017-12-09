#  Copyright (C) 2016-2017 Edgar Costa, Jeroen Sijsling
#  See LICENSE file for license details.

from sage.all import CRT_list, FiniteField, Integers, Matrix, NumberField, PolynomialRing, QQ
from sage.all import is_pseudoprime, kronecker_delta, prod, randint, set_random_seed 


def random_split_primes(field, primes = 1, bits = 62, seed = None, relative_ext = None):
    # input:
    # * Number Field
    # * number of split primes in the field
    # * number of bits desired in each prime
    # * random "seed"
    # * optionally we might also require that a certain polynomial splits completely over the split primes
    # Output a list [(p, root)] where:
    # * p is a partially split prime in the number field (and fully split in the relative extension)
    # * a root of the defining polynomial mod p
    lower_bound = 2 ** (bits - 1);
    if seed is None:
        k = randint(1,2 ** (bits - 2) )
    else:
        k = 0;
    p = lower_bound + 1 + 2*k;
    if field is QQ:
        f = [0, 1]
    else:
        f = field.defining_polynomial().list()

    output = [];
    ps = [];

    if relative_ext is not None:
        f_relative = relative_ext.list();


    while len(output) < primes:
        while not is_pseudoprime(p):
            p += 2;
        FF = FiniteField(p);
        FFz = PolynomialRing(FF, "z");
        roots = FFz(f).roots()
        if len(roots) > 0: #== f.degree():
            root = roots[0][0].lift();
            if relative_ext is not None:
                if len(FFz( reduce_list_split(f_relative, (p, root)) ).roots()) == len(f_relative) - 1:
                    output.append((p, root));
                    ps.append(p);
            else:
                output.append((p, root));
                ps.append(p);
        p += 2;
    return output;

def fractional_CRT_QQ(residues, ps_roots):
    residues_ZZ = [ r.lift() for r in residues];
    primes = [p for p, _ in ps_roots]
    lift = CRT_list(residues_ZZ, primes);
    N = prod(primes);
    M = Matrix([[1 ,lift, N]]);
    short_vector = Matrix(M.transpose().kernel().basis()).LLL()[0]
    return -short_vector[0]/short_vector[1];

def fractional_CRT_split(residues, ps_roots, K, BI_coordinates = None):
    if K is QQ:
        return fractional_CRT_QQ(residues, ps_roots);

    OK = K.ring_of_integers();
    BOK = OK.basis();

    residues_ZZ = [ r.lift() for r in residues];
    primes = [p for p, _ in ps_roots];
    lift = CRT_list(residues_ZZ, primes);
    lift_coordinates = OK.coordinates(lift);

    if BI_coordinates is None:
        p_ideals = [ OK.ideal(p, K.gen() - root) for p, root in ps_roots ];
        I = prod(p_ideals);
        BI = I.basis(); # basis as a ZZ-module
        BI_coordinates = [ OK.coordinates(b) for b in BI ];

    M = Matrix(Integers(), [ [ kronecker_delta(i,j) for j,_ in enumerate(BOK) ] + [ lift_coordinates[i] ] + [ b[i] for b in BI_coordinates ] for i in range(len(BOK)) ])
    # v = short_vector
    Kernel_Basis =  Matrix(Integers(), M.transpose().kernel().basis())
    v =Kernel_Basis.LLL()[0];
    #print v[:len(BOK)]
    if v[len(BOK)] == 0:
        return 0;
    return (-1/v[len(BOK)]) * sum( v[i] * b for i, b in enumerate(BOK));

def reduce_constant_split(x, p_root):
    x_list = list(x);
    p , root = p_root;
    return  FiniteField(p)(sum(xi * root**i for i, xi in enumerate(x_list)))


def reduce_constant_split_relative(x, p_root, relative_root):
    x_list = list(x);
    p , _ = p_root;
    return  FiniteField(p)(sum(reduce_constant_split(xi, p_root) * relative_root**i for i, xi in enumerate(x_list)))

def reduce_list_split( x_list, p_root):
    return [ reduce_constant_split(x, p_root) for x in x_list ]

def reduce_list_split_relative( x_list, p_root, relative_root):
    return [ reduce_constant_split_relative(x, p_root, relative_root) for x in x_list ]

def reduce_matrix_split(M, p_root):
    return Matrix([[ reduce_constant_split(x, p_root) for x in row ] for row in M.rows()]);

def reduce_matrix_split_relative(M, p_root, relative_root):
    return Matrix([[ reduce_constant_split_relative(x, p_root, relative_root) for x in row ] for row in M.rows()]);

def test_fractional_CRT(field, bits = 63, primes = 50):
    K = field;
    x = K.random_element();

    ps_roots = random_split_primes(field = K, bits = bits, primes = primes);

    residues = [ reduce_constant_split(x, p_root) for p_root in ps_roots];

    OK = K.ring_of_integers();
    p_ideals = [ OK.ideal(p, K.gen() - root) for p, root in ps_roots ];
    I = prod(p_ideals);
    if K is QQ:
        BI_coordinates = None;
    else:
        BI = I.basis(); # basis as a ZZ-module
        BI_coordinates = [ OK.coordinates(b) for b in BI ];
    x_recovered = fractional_CRT_split(residues, ps_roots, K, BI_coordinates);
    return x, x_recovered

def test_fractional_CRT_sample():
    QQt = PolynomialRing(QQ, "t");
    t = QQt.gen();
    for K in [QQ, NumberField(t**2  - 2, "a"),  NumberField(t**2  + 3, "a"),  NumberField(t**3  - 2, "a"),  NumberField(t**6 - t**5 + 2*t**4 + 8*t**3 - t**2 - 5*t + 7, "a"), NumberField(t**7 - 2*t**6 + 3*t**5 - 8*t**4 + 12*t**3 - 13*t**2 + 14*t - 6,"a"),NumberField(t**7 - 4*t**5 - t**4 - 5*t**3 + 4*t**2 - 4*t + 1,"a"), NumberField(t**5 - t**4 - 4*t**3 + 3*t**2 + 3*t - 1,"a")]:
        print "Testing fractional CRT on the %s" % K
        for i in range(100):
            set_random_seed(i)
            x, x_CRT = test_fractional_CRT(K)
            if x != x_CRT:
                print "Failed at i = %d" % i
                break;
        else:
            print "Test passed!\n"

