#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

import sys
from sage.all import CC, Infinity, Matrix, PolynomialRing, QQ
from sage.all import cputime, norm, vector, walltime
from misc import NF_embedding, toCCap, toCCap_list
from AbelJacobi import AJ1_digits

def compute_alpha_point(g, iaj, alpha, P0, digits, power, verbose, aggressive = True, append = ""):
    # input: 
    # * g, defining polynomial of the hyperelliptic curve
    # * iaj = InvertAJglobal class
    # * alpha = analytic representation of the endomorphism acting on the tangent space
    # * digits = how many digits to work with
    # * power, we will divide by 2**power before inverting the abel jacobi map
    # * P0, a point in the curve for which we will compute alpha*2(P - \infty) = P1 + P2 - \infty
    # * verbose
    # * agressive, raise ZeroDivisionError  if we get a P1[0] = P2[0] or one of the points is infty
    # * append, a string to append to the numberfield generators
    #
    # output: a dictionary with the following keys
    # TODO add the keys

    output = {};
    K = alpha.base_ring();
    CCap = iaj.C
    prec = CCap.precision()
    output['prec'] = prec;
    x0, y0 = P0;
    Kz = PolynomialRing(K, "z");
    z = Kz.gen();
    if verbose:
        print "compute_alpha_point()"
        print "P0 = %s" % (P0, )
    #deal with the field of definition
    if y0 in K:
        if K is QQ:
            L = K
            from_K_to_L = QQ.hom(1,QQ)
        else:
            L = K.change_names("b" + append)
            from_K_to_L = L.structure()[1];
    else:
        L, from_K_to_L = (z**2 - g(x0)).splitting_field("b" + append, simplify_all = True, map = True);
    
    output['L'] = L;
#    output['L_str'] = sage_str_numberfield(L, 'x','b'+append); 
    output['from_K_to_L'] = from_K_to_L
#    output['L_gen'] =  toCCap(L.gen(), 53);


    # figure out y0 in L
    y0_ap = toCCap(y0, prec);
    y0s = [elt for elt, _ in (z**2 - g(x0)).roots(L)];
    assert len(y0s) == 2;
    y0s_ap = [toCCap(elt, prec) for elt in y0s];
    if norm(y0_ap - y0s_ap[0]) <   norm(y0_ap - y0s_ap[1]):
        y0 = y0s[0];
    else:
        y0 = y0s[1]
    
    P0 = vector(L, [x0, y0])

    alpha_L = Matrix(L, [[from_K_to_L(elt) for elt in row] for row in alpha.rows()]);
    alpha_ap = Matrix(CCap, [toCCap_list(row, prec) for row in alpha_L.rows()]);
    
    output['P'] = P0;
    output['alpha'] = alpha_L; 
    

    if verbose:
        print "L = %s" % L
        print "%s = %s" % (L.gen(), toCCap(L.gen(), 53));
        print "P0 = %s" % (P0,)
        print "alpha = %s" % ([[elt for elt in row] for row in alpha_L.rows()], )


    P0_ap =  vector(CCap, toCCap_list(P0, prec + 192))

    if verbose:
        print "P0_ap = %s" % (vector(CC, P0_ap),)

    if verbose: 
        print "Computing AJ of P0..."
        c, w = cputime(), walltime()
        ajP0 = vector(CCap, AJ1_digits(g, P0_ap, digits));
        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
        print
    else:
        ajP0 = vector(CCap, AJ1_digits(g, P0_ap , digits))
     
    output['ajP'] = ajP0;

    aj2P0 = 2* ajP0;
    alpha2P0_an = alpha_ap * aj2P0
    if verbose:
        print "Working over %s" % CCap
    
    invertAJ_tries = 3;
    for i in range(invertAJ_tries):
        #try invertAJ_tries times to invertAJ
        try:
            if verbose:
                print "\n\ninverting AJ..."
                c, w = cputime(), walltime()
                alpha2P0_div = iaj.invertAJ(alpha2P0_an, power)
                print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
                print "\n\n"
            else:
                alpha2P0_div = iaj.invertAJ(alpha2P0_an, power)
            break;
        except AssertionError as error:
            print error
            if verbose:
                print "an assertion failed while inverting AJ.."
            if i == invertAJ_tries - 1:
                if verbose:
                    print "retrying with a new P0"
                    print;
                raise ZeroDivisionError;
            else:
                iaj.iajlocal.set_basepoints();
                if verbose:
                    print "retrying again with a new set of base points"

    
    if verbose:
        print "Computing the Mumford coordinates of alpha(2*P0 - 2*W)"
        c, w = cputime(), walltime()
        R0, R1 = alpha2P0_div.coordinates();
        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
    else:
        R0, R1 = alpha2P0_div.coordinates();


    points = [R0, R1]
    x_poly =  alpha2P0_div.x_coordinates()
    output['R'] = points;
    output['x_poly'] = x_poly;
    
    if (R0 in  [+Infinity, -Infinity]) or (R1 in  [+Infinity, -Infinity]):
        if aggressive:
            # we want to avoid these situations
            if verbose:
                print "One of the coordinates is at infinity"
                if R0 in  [+Infinity, -Infinity]:
                    print "R0 = %s" % (R0,)
                else:
                    print "R1 = %s" % (R1,)
                print "retrying with a new P0"
                print
            raise ZeroDivisionError;
        else:
            return output
    
    buffer = "# R0 = %s\n" % (vector(CC, R0),)
    buffer += "# R1 = %s\n" % (vector(CC, R1),)
    buffer += "# and \n# x_poly = %s\n" % (vector(CC, x_poly),)
    
    if verbose:
        print buffer
    assert len(x_poly) == 3


    algx_poly = [NF_embedding(coeff, L) for coeff in x_poly]
    buffer = "algx_poly = %s;\n" % (algx_poly,)
    if L != QQ:
        buffer += "#where %s ~ %s\n" % (L.gen(), L.gen().complex_embedding())
    buffer += "\n" 
    if verbose:
        print buffer
        sys.stdout.flush()
        sys.stderr.flush()
    
    output['algx_poly'] = algx_poly;

    if  None in  algx_poly:
        if aggressive:
            if verbose:
                print "No algebraic expression  for the polynomial"
                print "retrying with a new P0"
                sys.stdout.flush()
                sys.stderr.flush()
            raise ZeroDivisionError;
        else:
            return output;
    #c, b, a = algx_poly
    #if aggressive and b**2 - 4*a*c == 0:
    #    raise ZeroDivisionError
    if verbose:
        print "Done compute_alpha_point()"
    return output

# nonsense from the past
#    if False:
#        # some nonsense from the past         
#        M = Lw(algx_poly).splitting_field('c'+append, simplify_all = True);
#        
#        output['M'] = M;
#        if verbose:
#            print "M is the field where the x coordinates of R0 and R1 are defined"
#            print "M = %s" % M 
#
#        alg_points = [None, None];
#        for i, P in enumerate(points):
#            alg_x = NF_embedding(P[0], M);
#            if alg_x is None or g(alg_x) == 0:
#                if aggressive:
#                    raise ZeroDivisionError
#                else:
#                    return output
#            P = [ toCCap(alg_x, prec), P[1] ];
#            
#            P_branch = branch(g, P);
#            alg_points[i] = (alg_x, P_branch)
#            buffer += "algR%d_x = %s;\n" % (i, alg_x);
#            buffer += "R%d_branch = %s;\n" % (i, P_branch);
#            buffer += "\n"
#        if alg_points[0][0] == alg_points[1][0] and aggressive:
#            raise ZeroDivisionError
#
#        output['algR'] = alg_points;
#        if M is not QQ:
#            buffer += "#where %s ~ %s" % (M.gen(), M.gen().complex_embedding())
#        if verbose:
#            print buffer


