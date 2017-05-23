#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

from bound_rosati import bound_rosati
from sage.all import CC, ComplexField, Matrix, NumberField, PolynomialRing,  QQ 
from sage.all import copy, ceil, norm, sqrt, prod, vector 
from misc import almost_equal,  not_equal,  toCCap
from newton_lift import  trace_and_norm_ladic
from verify_algebraically import verify_algebraically

def add_trace_and_norm_ladic(g, D, alpha_geo, verbose = True):
    if verbose:
        print "add_trace_and_norm_ladic()"
    #load Fields
    L = D['L']
    if L is QQ:
        QQx = PolynomialRing(QQ, "x");
        L = NumberField(QQx.gen(), "b")
    
    prec = D['prec']
    CCap = ComplexField(prec);

    # load endo
    alpha = D['alpha']
    rosati = bound_rosati(alpha_geo)

    if alpha.base_ring() is not L:
        alpha_K = copy(alpha)
        alpha = Matrix(L, 2,2)
        #shift alpha field from K to L
        for i, row in enumerate(alpha_K.rows()):
            for j, elt in enumerate(row):
                if elt not in L:
                    assert elt.base_ring().absolute_polynomial() == L.absolute_polynomial()
                    alpha[i, j] = L(elt.list())
                else:
                    alpha[i, j] = L(elt)



    # load algx_poly
    algx_poly_coeff = D['algx_poly']


    #sometimes, by mistake the algx_poly is defined over K where K == L, but with a different name
    for i, elt in enumerate(algx_poly_coeff):
        if elt not in L:
            assert elt.base_ring().absolute_polynomial() == L.absolute_polynomial()
            algx_poly_coeff[i] = L(elt.list())
        else:
            algx_poly_coeff[i] = L(elt)

    x_poly = vector(CCap, D['x_poly']);
    for i in [0,1]:
        assert almost_equal(x_poly[i] , algx_poly_coeff[i]), "%s != %s" % (algx_poly_coeff[i], x_poly[i])

    # load P
    P0 = vector(L, [D['P'][0], D['P'][1]])
    for i, elt in enumerate(P0):
        if elt not in L:
            assert elt.base_ring().absolute_polynomial() == L.absolute_polynomial()
            P0[i] = L(elt.list())
        else:
            P0[i] = L(elt)
    if verbose:
        print "P0 = %s" % (P0,);

    # load image points, P1 and P2

    L_poly = PolynomialRing(L, "xL");
    xL = L_poly.gen();

    Xpoly = L_poly(algx_poly_coeff)
    if Xpoly.is_irreducible():
        M = Xpoly.root_field("c")
    else:
        # this avoids bifurcation later on in the code, we don't want to be always checking if M is L
        M = NumberField(xL, "c")
    
    # trying to be sure that we keep the same complex_embedding...
    M_complex_embedding = 0 
    if L.gen() not in QQ:
        M_complex_embedding = None
        Lgen_CC = toCCap(L.gen(), prec = prec)
        for i, _ in enumerate(M.complex_embeddings()):
            if norm(Lgen_CC-  M(L.gen()).complex_embedding(prec = prec, i = i)) < CCap(2)**(-0.7 * Lgen_CC.prec()) * Lgen_CC.abs():
                M_complex_embedding = i;

        assert M_complex_embedding is not None, "\nL = %s\n%s = %s\n%s" % (L, L.gen(), Lgen_CC, M.complex_embeddings())
    
    M_poly = PolynomialRing(M, "xM");
    xM = M_poly.gen();

    # convert everything to M
    P0_M = vector(M, [elt for elt in P0]);
    alpha_M = Matrix(M, [[elt for elt in row] for row in alpha.rows()])
    Xpoly_M = M_poly(Xpoly) 

    for i in [0,1]:
        assert almost_equal(x_poly[i] , Xpoly_M.list()[i], ithcomplex_embedding = M_complex_embedding), "%s != %s" % (Xpoly_M.list()[i], x_poly[i])

    P1 =  vector(M, 2)
    P1_ap = vector(CCap, D['R'][0])
    P2 =  vector(M, 2)
    P2_ap = vector(CCap, D['R'][1])



    M_Xpoly_roots = Xpoly_M.roots()
    assert len(M_Xpoly_roots) > 0

    Ypoly_M = prod([xM**2 - g(root) for root,_ in M_Xpoly_roots]) # also \in L_poly
    


    assert sum(m for _, m in Ypoly_M.roots(M)) == Ypoly_M.degree(), "%s != %s\n%s\n%s" % (sum(m for _, m in Ypoly_M.roots(M)), Ypoly_M.degree(), Ypoly_M, Ypoly_M.roots(M))



    if len(M_Xpoly_roots) == 1:
        # we have a double root
        P1[0] = M_Xpoly_roots[0][0]
        P2[0] = M_Xpoly_roots[0][0]
        ae_prec = prec*0.4
    else:
        assert len(M_Xpoly_roots) == 2
        ae_prec = prec
        # we have two distinct roots
        P1[0] = M_Xpoly_roots[0][0]
        P2[0] = M_Xpoly_roots[1][0]
        if not_equal(P1_ap[0], P1[0], ithcomplex_embedding = M_complex_embedding):
            P1[0] = M_Xpoly_roots[1][0]
            P2[0] = M_Xpoly_roots[0][0]
    
    
    assert almost_equal(P1_ap[0], P1[0], ithcomplex_embedding = M_complex_embedding, prec = ae_prec), "\n%s = %s \n != %s" % (P1[0], toCCap(P1[0], ithcomplex_embedding = M_complex_embedding), CC(P1_ap[0]))
    assert almost_equal(P2_ap[0], P2[0], ithcomplex_embedding = M_complex_embedding, prec = ae_prec), "\n%s = %s \n != %s" % (P2[0], toCCap(P2[0], ithcomplex_embedding = M_complex_embedding), CC(P2_ap[0]))

    # figure out the right square root

    # pick the default branch
    P1[1] = sqrt(g(P1[0]))
    P2[1] = sqrt(g(P2[0]))
    
    if 0 in [P1[1], P2[1]]:
        print "one of image points is a Weirstrass point"
        print P1
        print P2
        raise ZeroDivisionError

    #switch if necessary
    if not_equal(P1_ap[1], P1[1], ithcomplex_embedding = M_complex_embedding, prec = ae_prec):
        P1[1] *= -1

    if not_equal(P2_ap[1], P2[1], ithcomplex_embedding = M_complex_embedding, prec = ae_prec):
        P2[1] *= -1

    # double check
    for i in [0,1]:
        assert almost_equal(P1_ap[i] , P1[i], ithcomplex_embedding = M_complex_embedding, prec = ae_prec), "%s != %s" % (P1_ap[i], P1[i])
        assert almost_equal(P2_ap[i] , P2[i], ithcomplex_embedding = M_complex_embedding, prec = ae_prec), "%s != %s" % (P2_ap[i], P2[i])

    # now alpha, P0 \in L
    # P1, P2 \in L

    if verbose:
        print "P1 = %s\nP2 = %s" % (P1, P2)
        print "Computing the trace and the norm ladically\n"
        trace_and_norm = trace_and_norm_ladic(L, M, P0_M, P1, P2, g, 2*alpha_M, 16*rosati, primes = ceil(prec * L.degree()/61));
    else:
        trace_and_norm = trace_and_norm_ladic(L, M, P0_M, P1, P2, g, 2*alpha_M, 16*rosati, primes = ceil(prec * L.degree()/61))

    # Convert the coefficients to polynomials
    trace_numerator, trace_denominator, norm_numerator, norm_denominator = [L_poly(coeff) for coeff in trace_and_norm];

    assert trace_numerator(P0[0]) == trace_denominator(P0[0]) * -algx_poly_coeff[1], "%s/%s (%s) != %s" % (trace_numerator, trace_denominator, P0[0], -algx_poly_coeff[1])
    assert norm_numerator(P0[0]) == norm_denominator(P0[0]) * algx_poly_coeff[0], "%s/%s (%s) != %s" % (norm_numerator, norm_denominator, P0[0], algx_poly_coeff[0])
    buffer =  "# x1 + x2 = degree %d/ degree %d\n" % (trace_numerator.degree(), trace_denominator.degree())
    buffer += "# = (%s) / (%s) \n" %  (trace_numerator, trace_denominator)
    buffer += "# max(%d, %d) <= %d\n\n" % (trace_numerator.degree(), trace_denominator.degree(), 16 * rosati);

    buffer += "# x1 * x2 = degree %d/ degree %d\n" % (norm_numerator.degree(), norm_denominator.degree())
    buffer += "# = (%s) / (%s) \n" %  (norm_numerator, norm_denominator)
    buffer += "# max(%d, %d) <= %d\n" % (norm_numerator.degree(), norm_denominator.degree(), 16 * rosati);
    

    if verbose:
        print buffer;
        print "\n"
    assert max(trace_numerator.degree(), trace_denominator.degree()) <= 16 * rosati
    assert max(norm_numerator.degree(), norm_denominator.degree()) <= 16 * rosati

    if verbose:
        print "Veritfying if x1*x2 and x1 + x2 are correct..."

    verified = verify_algebraically(g, P0, alpha, trace_and_norm, verbose = verbose)
    if verbose:
        print "\nDoes it act on the tangent space as expected? %s\n" % verified
        print "Done add_trace_and_norm_ladic()"
 
    return verified, [trace_numerator.list(), trace_denominator.list(), norm_numerator.list(), norm_denominator.list()]



# DEPRECATED
# and somewhat broken
#def add_trace_and_norm_ap(D, verbose = True):
#    for t in range(10):
#        try:
#            alpha = D['alpha']
#            alpha_geo = D['alpha_geo']
#
#
#            # 4 from the doubling
#            # 2 for the PS terms
#            rosati = bound_rosati(alpha_geo)
#            if verbose:
#                print "alpha_geo =\n%s "  % alpha_geo
#            total_degree= 16  * rosati;
#            hard_bound =  2 * total_degree + 1;
#            soft_bound = hard_bound + 10;
#            working_precision = max(soft_bound * 10 * 2**(t + 1), 2000);
#            CCap = ComplexField(working_precision);
#            buffer = "# bound_rosati(alpha_geo) = %d\n" % rosati;
#            buffer += "# hard_bound = 2 * 16 * bound_rosati(alpha_geo) + 1 = %d;\n" % hard_bound
#            #buffer = "data[%d]['hard_bound'] = %d\n;" % (t, hard_bound)
#            buffer += "# soft_bound = hard_bound + 10 ;\n"
#            #buffer = "data[%d]['soft_bound'] = %d\n;" % (t, soft_bound)
#            buffer += "data[%d]['lift_working_precision'] = %d;\n" % (i, working_precision, )
#            if verbose:
#                print buffer
#
#            algR0, algR1 = D['algR'];
#            R0, R1 = D['R'];
#            algx_poly = D['algx_poly'];
#
#
#            L = D['L']; # the numberfield where P and algx_poly are defined
#            LW = PowerSeriesRing(L, "W", default_prec = soft_bound);
#            W = LW.gen();
#            L = PolynomialRing(L, "T");
#            T = L.gen();
#
#
#            PS = PowerSeriesRing(CCap, "Z", default_prec = soft_bound)
#            Z = PS.gen();
#            #P.<U> = CCap[]
#
#            P0 = vector(L, D['P'])
#            P0_an = vector(CCap, [toCCap(x, working_precision) for x in P0])
#
#
#
#            x1_an = toCCap(algR0[0], working_precision);
#            b1 = branch(g, vector(CCap, [x1_an, R0[1]]));
#            P1_an = vector(CCap, [x1_an, b1*sqrt(g(x1_an))])
#
#
#            x2_an = toCCap(algR1[0], working_precision)
#            b2 = branch(g, vector(CCap, [x2_an, R1[1]]));
#            P2_an = vector(CCap, [x2_an, b2*sqrt(g(x2_an))])
#
#            alpha_ap = Matrix(CCap, [[ toCCap(elt, working_precision) for elt in row] for row in alpha.rows()])
#
#            if verbose:
#                print "lifting over CC"
#                c, w = cputime(), walltime()
#                x1, x2 = newton_linear(P0_an, P1_an, P2_an, g, 2*alpha_ap, PS, soft_bound)
#                print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
#            else:
#                x1, x2 = newton_linear(P0_an, P1_an, P2_an, g, 2*alpha_ap, PS, soft_bound)
#            trace_approx = x1 + x2
#            norm_approx = x1 * x2
#            trace_series = NF_series_embedding(trace_approx, LW);
#            norm_series = NF_series_embedding(norm_approx, LW);
#
#
#            trace_numerator, trace_denominator = rational_reconstruct_poly(LT(trace_series), T**trace_series.prec())
#            if trace_numerator.degree() == (trace_series.prec() - 1)//2 or trace_denominator.degree() == (trace_series.prec() - 1)//2:
#                print "Error: rational reconstruction might have failed"
#                print "trace_series.prec() = %d\ntrace_numerator.degree() = %d\ntrace_denominator.degree() = %d" % (trace_series.prec(), trace_numerator.degree(), trace_denominator.degree())
#                raise ZeroDivisionError;
#
#            norm_numerator, norm_denominator = rational_reconstruct_poly(LT(norm_series), T**norm_series.prec())
#            if norm_numerator.degree() == (norm_series.prec() - 1)//2 or norm_denominator.degree() == (norm_series.prec() - 1)//2:
#                print "Error: rational reconstruction might have failed"
#                print "norm_series.prec() = %d\norm_numerator.degree() = %d\norm_denominator.degree() = %d" % (norm_series.prec(), norm_numerator.degree(), norm_denominator.degree())
#                raise ZeroDivisionError;
#
#
#            trace_Z = PS( toCCap_list( (LW(trace_numerator)/LW(trace_denominator)).list(), working_precision) )
#            norm_Z = PS( toCCap_list( (LW(norm_numerator)/LW(norm_denominator)).list(), working_precision) )
#            trace_diff = (trace_Z - trace_approx).padded_list();
#            norm_diff =  (norm_Z - norm_approx).padded_list();
#            trace_linf = max([CC(elt.abs()) for elt in trace_diff ])
#            norm_linf = max([CC(elt.abs()) for elt in norm_diff ])
#            # Utrace_numerator, Utrace_denominator = (x1*x2).pade(total_degree, total_degree)
#            # print vector(CC,Utrace_numerator.list())
#            # print vector(CC,Utrace_denominator.list())
#
#
#            # the series was based at P0[0], now we can change it to 0
#            trace_numerator = trace_numerator.subs(T = T - P0[0])
#            trace_denominator = trace_denominator.subs(T = T - P0[0])
#            norm_numerator = norm_numerator.subs(T = T - P0[0])
#            norm_denominator = norm_denominator.subs(T = T - P0[0])
#
#
#            assert trace_numerator(P0[0]) == trace_denominator(P0[0]) * -algx_poly[1], "%s/%s (%s) != %s" % (trace_numerator, trace_denominator, P0[0], -algx_poly[1])
#            assert norm_numerator(P0[0]) == norm_denominator(P0[0]) * algx_poly[0], "%s/%s (%s) != %s" % (norm_numerator, norm_denominator, P0[0], algx_poly[0])
#
#
#            buffer = "# x1 + x2 = degree %d/ degree %d\n" % (trace_numerator.degree(), trace_denominator.degree())
#            buffer += "# = (%s) / (%s) \n" %  (trace_numerator, trace_denominator)
#            buffer += "# %d + %d <= %d\n\n" % (trace_numerator.degree(), trace_denominator.degree(), total_degree);
#            buffer += "# | x1 + x2 - (x1 + x2 + O(T**%d)  )_approx|_inf = \ndata[%d]['trace_linf'] = %s\n" % ( trace_series.prec(), i, RR( trace_linf ));
#            buffer += "# x1 * x2 = degree %d/ degree %d\n" % (norm_numerator.degree(), norm_denominator.degree())
#            buffer += "# = (%s) / (%s) \n" %  (norm_numerator, norm_denominator)
#            buffer += "# %d + %d <= %d\n" % (norm_numerator.degree(), norm_denominator.degree(), total_degree);
#            buffer += "# | x1 * x2 - (x1 * x2 + O(T**%d)  )_approx|_inf = \ndata[%d]['norm_linf'] = %s\n" % ( norm_series.prec(), i, RR( norm_linf ));
#            buffer += "# x1 + x2 = \n";
#            buffer += "# %s;\n" % ( trace_numerator.list(),);
#            buffer += "# %s;\n" % ( trace_denominator.list(),);
#            buffer += "# x1 * x2 = \n";
#            buffer += "# %s;\n" % ( norm_numerator.list(),);
#            buffer += "# %s;\n" % ( trace_denominator.list(),);
#           
#            assert trace_numerator.degree() + trace_denominator.degree() <= 16 * rosati
#            assert norm_numerator.degree() + norm_denominator.degree() <= 16 * rosati
#            if verbose:
#                print buffer;
#
#
#            verified = verify_algebraically(g, P0, alpha,  trace_and_norm)
#            if verbose:
#                print "# Does it act on the tangent space as expected? %s" % verified
#            
# 
#            break
#
#        except ZeroDivisionError as err:
#            print("ZeroDivisionError: {0}".format(err))
#            print "Not working with enough precision";
#            pass;
#        except ValueError as err:
#            print("ValueError: {0}".format(err))
#            print "Not working with enough precision?";
#            pass;
#    else:
#        print "something went wrong..."
#        raise ZeroDivisionError
#    return  verified, [elt.list() for elt in [trace_numerator, trace_denominator, norm_numerator, norm_denominator]]
