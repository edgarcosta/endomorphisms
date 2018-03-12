#  Copyright (C) 2016, Edgar Costa
#  See LICENSE file for license details

from sage.all import FiniteField, O, PolynomialRing, PowerSeriesRing, QQ, ZZ
from sage.all import cputime, norm, prod, sqrt, vector, walltime
from FractionalCRT import fractional_CRT_split, random_split_primes,reduce_constant_split, reduce_list_split, reduce_list_split_relative, reduce_matrix_split_relative
from rational_reconstruct import rational_reconstruct_poly

def trace_and_norm_ladic(L, M, P0, P1, P2, f, alpha, degree_bound, primes = 120, bits = 62, verbose = True):
    if verbose:
        print "trace_and_norm_ladic()"
        print "primes = %d, bits = %d" % (primes, bits)
    # Input:
    # * L the number field where P0, alpha, norm and trace are defined over
    # * M a relative extension of L where P1 and P2 are defined over
    # * P0 initial point already embedded in M
    # * P1, P2 image points alpha(2 * P0 - \infty) = P1 + P2 - \infty
    # * alpha the matrix represing the endomorphism on the tangent space in M
    # * f the defining polynomial of the curve such that y^2 = f(x)
    # * degree bound for the trace and norm as rational maps
    # Note: Due to inconsistencies of the complex embeddings of L and M we require P0 and alpha to be given in M
    # Output:
    # * [trace_numberator, trace_denominator, norm_numberator, norm_denominator]
    # represented as lists of elements in L

    
    for arg in [P0, P1, P2, alpha]:
        assert arg.base_ring() is M

    hard_bound = 2*degree_bound + 1
    soft_bound = hard_bound + 10;
    if verbose:
        print "hard_bound = %d\nsoft_bound = %d" % (hard_bound,soft_bound)
        print "Generating split primes..."
        c, w = cputime(), walltime()
        ps_roots = random_split_primes(field = L, bits = bits, primes = primes, relative_ext = M.relative_polynomial());
        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
    else:
         ps_roots = random_split_primes(field = L, bits = bits, primes = primes, relative_ext = M.relative_polynomial());




    to_lift_raw = [];
    if verbose:
        print "Applying newton lift at each prime"
    for p_root in ps_roots:
        ell, root = p_root;

        FF = FiniteField( ell );
        FF_xell = PolynomialRing(FF, "xell");

        Lpoly = FF_xell(reduce_list_split(M.relative_polynomial().list(), p_root));
        assert Lpoly.degree() == len(Lpoly.roots());
        Lroot = Lpoly.roots()[0][0].lift();

        PSring_ell = PowerSeriesRing(FF, "Tell", default_prec = soft_bound )

        P0_ell = vector(FF, reduce_list_split_relative(P0, p_root, Lroot));
        P1_ell = vector(FF, reduce_list_split_relative(P1, p_root, Lroot));
        P2_ell = vector(FF, reduce_list_split_relative(P2, p_root, Lroot));

        f_ell = FF_xell(reduce_list_split(f.list(), p_root));
        alpha_ell = reduce_matrix_split_relative(alpha, p_root, Lroot);


        x1_ell, x2_ell = newton_linear(P0_ell, P1_ell, P2_ell, f_ell, alpha_ell, PSring_ell, soft_bound);

        trace_series_ell = x1_ell + x2_ell

        norm_series_ell = x1_ell * x2_ell

        trace_numerator_ell, trace_denominator_ell = rational_reconstruct_poly(FF_xell(trace_series_ell.list()), FF_xell.gen()**trace_series_ell.prec())
        norm_numerator_ell, norm_denominator_ell = rational_reconstruct_poly(FF_xell(norm_series_ell.list()), FF_xell.gen()**norm_series_ell.prec())

        # shift the series back to zero
        shift = FF_xell.gen() - P0_ell[0]

        factors_ell = [trace_numerator_ell(shift).list(), trace_denominator_ell(shift).list(), norm_numerator_ell(shift).list(), norm_denominator_ell(shift).list()];
        degrees_ell = tuple([ len(elt) for elt in factors_ell ]);
        to_lift_raw.append((p_root, degrees_ell, factors_ell))
    if verbose:
        print "Getting rid of bad primes..."
    # get rid of badprimes
    degrees_count = {};
    for _, degrees_ell, _ in to_lift_raw:
        if degrees_ell in degrees_count:
            degrees_count[degrees_ell] += 1
        else:
            degrees_count[degrees_ell] = 1
    max_value = 0;
    max_arg = None;

    for degrees_ell, count in degrees_count.iteritems():
        if count > max_value:
            max_arg = degrees_ell;
            max_value = count;

    output = [None] * 4;
    for i, elt in enumerate(max_arg):
        output[i] = [0]*elt

    to_lift = [];
    ps_roots = [];

    for p_root, degrees_ell, factors_ell in to_lift_raw:
        if degrees_ell == max_arg:
            ps_roots.append(p_root);
            to_lift.append(factors_ell);

    extra_factor = to_lift[-1];
    extra_p_root = ps_roots[-1];


    to_lift = to_lift[:-1];
    ps_roots = ps_roots[:-1];

    OL = L.ring_of_integers();
    p_ideals = [ OL.ideal(p, L.gen() - r) for p, r in ps_roots ];
    I = prod(p_ideals);
    if L is QQ:
        BI_coordinates = None;
    else:
        BI = I.basis(); # basis as a ZZ-module
        BI_coordinates = [ OL.coordinates(b) for b in BI ];
    
    if verbose:
        print "Lifting everything"
    for i, lift in enumerate(output):
        for j, elt in enumerate(lift):
            residues = [ residue[i][j] for residue in to_lift ]
            output[i][j] = fractional_CRT_split(residues = residues, ps_roots = ps_roots, K = L, BI_coordinates = BI_coordinates)
            assert reduce_constant_split( output[i][j], extra_p_root) == extra_factor[i][j]
    if verbose:
        print "trace_and_norm_ladic() Done"

    return output



def newton_linear(P0, P1, P2,  f, alpha, PSring, prec):
    assert P0.base_ring() is P1.base_ring()
    assert P1.base_ring() is P2.base_ring()
    assert P0.base_ring() is alpha.base_ring()
    assert P0.base_ring() is PSring.base_ring()
    assert f.base_ring() in [QQ, ZZ, P0.base_ring()]
    T = PSring.gen();
    #xi_0 = xi(0)
    x0_0 = PSring(P0[0]);
    x1_0 = PSring(P1[0]);
    x2_0 = PSring(P2[0]);


    y0_0 = P0[1];
    y1_0 = P1[1];
    y2_0 = P2[1];

    b0, b1, b2 = 1, 1, 1;

    py0_0 = sqrt(f(x0_0 + T)).list()[0];
    py1_0 = sqrt(f(x1_0 + T)).list()[0];
    py2_0 = sqrt(f(x2_0 + T)).list()[0];


    if PSring.base_ring().is_exact():
        if y0_0 != py0_0:
            b0 = -1
        if y1_0 != py1_0:
            b1 = -1
        if y2_0 != py2_0:
            b2 = -1
        assert y0_0 == b0 * py0_0, "wrong branch at P0? %s != %s" % (y0_0, b0*py0_0);
        assert y1_0 == b1 * py1_0, "wrong branch at P1? %s != %s" % (y1_0, b1*py1_0);
        assert y2_0 == b2 * py2_0, "wrong branch at P2? %s != %s" % (y2_0, b2*py2_0);
    else:
        if norm(y0_0 - py0_0) > norm(y0_0 + py0_0):
            b0 = -1;
        if norm(y1_0 - py1_0) > norm(y1_0 + py1_0):
            b1 = -1;
        if norm(y2_0 - py2_0) > norm(y2_0 + py2_0):
            b2 = -1;
        closetozero = 2**(-0.8 * PSring.base_ring().prec())
        assert norm(y0_0 + b0*py0_0) > norm(y0_0), "%.3e vs %.3e wrong branch?" % (norm(y0_0 - b0*py0_0), norm(y0_0 + b0*py0_0) );
        assert norm(y0_0 - b0*py0_0) < norm(y0_0)*closetozero, "%.3e vs %.3e wrong branch?" % (norm(y0_0 - b0*py0_0), norm(y0_0 + b0*py0_0) );

        assert norm(y1_0 + b1*py1_0) > norm(y1_0), "%.3e vs %.3e wrong branch?" % (norm(y1_0 - b1*py1_0), norm(y1_0 + b1*py1_0));
        assert norm(y1_0 - b1*py1_0) < norm(y1_0)*closetozero, "%.3e vs %.3e wrong branch?" % (norm(y1_0 - b1*py1_0), norm(y1_0 + b1*py1_0));

        assert norm(y2_0 + b2*py2_0) > norm(y2_0), "%.3e vs %.3e wrong branch?" % (norm(y2_0 - b2*py2_0), norm(y2_0 + b2*py2_0) );
        assert norm(y2_0 - b2*py2_0) < norm(y2_0)*closetozero, "%.3e vs %.3e wrong branch?" % (norm(y2_0 - b2*py2_0), norm(y2_0 + b2*py2_0) );

    
    x1 = x1_0 + O(T);
    x2 = x2_0 + O(T);

    p = 1;

    Tsub = x0_0 + T;
    fps = PSring(f)(Tsub);
    sqrtfps = sqrt(fps);


    row0 = PSring(alpha.row(0).list())(Tsub);
    row1 = PSring(alpha.row(1).list())(Tsub);
    
    if x1_0 != x2_0:
    # Iteration to solve the differential equation, gains 1 term
    # alternatively we could solve a linear ODE to double the number of terms
        while p < prec:
            m = b0 / ( (x2 - x1)  * sqrtfps ) ;
            dx1 = m * b1 * sqrt(f(x1)) * ( row0 * x2 - row1 )
            dx2 = m * b2 * sqrt(f(x2)) * ( row1 - row0 * x1 )
            p += 1
            x1 = x1_0 + dx1.integral() + O(T**p)
            x2 = x2_0 + dx2.integral() + O(T**p)
    else:
        #the degenerate case
        row0_y = row0 * b0 / sqrtfps
        half = P0.base_ring()(1/2); 
        while p < prec:
            dx1 = half *  b1 * sqrt(f(x1)) * row0_y
            x1 = x1_0 + dx1.integral() + O(T**p)
            p += 1
        # check second equation
        if  2 * x1 * x1.derivative() * b1 / sqrt(f(x1)) != row1 *  b0 / sqrtfps:
            print "x1(0) = x2(0), but x1 != x2"
            raise ZeroDivisionError
        
        x2 = x1;



    return x1, x2

