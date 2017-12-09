#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

import sys, os
from sage.all import Matrix, PolynomialRing, QQ
from sage.all import ceil, cputime,  log, set_random_seed, walltime
import lmfdb
from misc import sage_str_numberfield, convert_magma_matrix_to_sage, find_rational_point
from InvertAJglobal import InvertAJglobal
from alpha_point import compute_alpha_point
from trace_and_norm import add_trace_and_norm_ladic


def process_curve_batch(folder, batch, digits = 600, power = 15):
#    cursor = list(C.genus2_curves.curves.find({'is_simple_geom': true,'real_geom_end_alg': u'R x R'}).sort([("cond", ASCENDING)]))
    totalcurves = len(batch)
    D = {};
    success = 0;
    total = 0;
    failed = [];
    for i, curve in enumerate(batch):
        label = curve['label']
        stdout_filename = os.path.join(folder, label + ".stdout")
        output_filename = os.path.join(folder, label + ".sage")

        print "%d of %d : label = %s" % (i+1, totalcurves, label)

        sys.stdout = open(stdout_filename, 'w')

        c, w = cputime(), walltime()
        data, sout = process_curve_lmfdb(curve, digits, power = power, verbose = True,  internalverbose = False)
        print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
        output_file = open(output_filename, 'w');
        output_file.write(sout);
        output_file.close()
        D[label] = True

        sys.stdout.flush()
        os.fsync(sys.stdout.fileno())
        sys.stdout = sys.__stdout__
        total += 1
        if label in D.keys():
            success += 1
            print "Done: label = %s" % (label)
            print os.popen("tail -n 1 %s" % stdout_filename).read()
        else:
            failed += [label];
            print "ERROR: label = %s\n" % (label)
        print "Success rate: %.0f%%\n" % (100.*success/total)
        print "Failed so far %s" %  failed

def process_curve_standalone(label, digits = 600,  power = 15, verbose = True, internalverbose = False, folder = None):
    import os
    set_random_seed(1);
    import random;
    random.seed(1)
    C = lmfdb.getDBconnection()
    curve = C.genus2_curves.curves.find_one({'label' : label})
    if curve is  None:
        print "Wrong label"
        return 2;
    endo_alg = curve['real_geom_end_alg'];
    # deals with the default folders
    if folder is None:
        if not '__datadir__' in globals():
            import os
            import inspect
            filename = inspect.getframeinfo(inspect.currentframe())[0];
            __datadir__ = os.path.dirname(filename) + "/data/" 
        base_folder = __datadir__
        if endo_alg == 'R':
            print "End = QQ, Nothing to do"
            return 0;
        elif endo_alg == 'R x R':
                if curve['is_simple_geom']:
                    type_folder = 'simple/RM';
                else:
                    type_folder = 'split/nonCM_times_nonCM' ;
        elif endo_alg == 'C x R':
            type_folder = 'split/CM_times_nonCM';
        elif endo_alg == 'C x C':
            if curve['is_simple_geom']:
                type_folder = 'simple/CM';
            else:
                type_folder = 'split/CM_times_CM';
        elif endo_alg == 'M_2(R)':
            if curve['is_simple_geom']:
                type_folder = 'simple/QM';
            else:
                type_folder = 'split/nonCM_square'
        elif endo_alg == 'M_2(C)':
            type_folder = 'split/CM_square'
        else:
            print "did I forget a case?"
            return 1;
        folder = os.path.join(base_folder, type_folder)
        if not os.path.exists(folder):
            print "Creating dir: %s" % folder
            os.makedirs(folder);

    assert os.path.exists(folder);
    #filenames
    stdout_filename = os.path.join(folder, label + ".stdout")
    stderr_filename = os.path.join(folder, label + ".stderr")
    output_filename = os.path.join(folder, label + ".sage")

    sys.stdout = open(stdout_filename, 'w', int(0))
    sys.stderr = open(stderr_filename, 'w')

    c, w = cputime(), walltime()
    data, sout, verified = process_curve_lmfdb(curve, digits, power = power, verbose = True,  internalverbose = False)
    print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
    output_file = open(output_filename, 'w');
    output_file.write(sout);
    output_file.close()

    sys.stdout.flush()
    os.fsync(sys.stdout.fileno())
    sys.stderr.flush()
    os.fsync(sys.stderr.fileno())
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print "Done: label = %s, verified = %s" % (label, verified)
    print os.popen("tail -n 1 %s" % stdout_filename).read()
    if verified:
        return 0;
    else:
        return 1;






    

def process_curve_lmfdb(curve, digits = 600, power = 15, verbose = True, internalverbose = False):
    #load(__endodir__ + 'Initialize.sage')
    label = curve['label']
    
    # Pick the model
    # minimal model
    fmin, hmin = eval(curve['eqn'])

    Qx = PolynomialRing(QQ, "x");
    g = 4*Qx(fmin) + Qx(hmin)**2
    return certify_heuristic(g, label, digits, power, verbose, internalverbose)

def certify_heuristic(g, label = None, digits = 600, power = 15, verbose = True, internalverbose = False):
    from heuristic_endomorphisms import EndomorphismData
    out = "";
    curve_dict = {};
    curve_dict['label'] = label
    curve_dict['digits'] = digits
    prec = ceil(log(10)/log(2)*digits)
    curve_dict['prec'] = prec;
    buffer = "curve_dict = {};\n";
    for key in ['digits', 'prec']:
        buffer += "curve_dict['%s'] = %s;\n" % (key, curve_dict[key]);
    buffer += "curve_dict['%s'] = '%s';\n" % ('label',label)
    if verbose:
        print buffer
    out += buffer;

    #Ambient ring
    buffer = "";
    buffer += "\n# basic rings\n";
    buffer += "curve_dict['%s'] = %s;\n" % ('QQx', 'PolynomialRing(QQ, \'x\')');
    buffer += "# curve_dict, x, and the number field's generator are the only variables are the only global variables\n\n"
    buffer += "x = curve_dict['QQx'].gen();\n\n";
    out+= buffer;

    Qx = PolynomialRing(QQ, "x")
    x = Qx.gen()
    curve_dict['QQx'] = Qx;
   
    # force g in Qx
    g = Qx(g.list())


    # Compute EndomorphismData
    # field of definition
    # and endomorphisms
    xsubs_list = [x, x + 1, x - 1, 1 - x, -1 - x, x + 2, x - 2, -x - 2, 2 - x]
    if label == '540800.a.540800.1':
        xsubs_list = xsubs_list[1:]
    for xsubs in xsubs_list:
        try:
            gtry = g(xsubs);
            if gtry.degree() == 5:
                gtry = Qx( gtry(x**(-1)) * x**6 )
            if gtry.degree() == 5:
                gtry = Qx( gtry((x - 1)/x) * x**6 )

            if verbose:
                print "Computing the EndomorphismData..."
                c, w = cputime(), walltime()
                End = EndomorphismData(gtry, prec = digits)
                print "Time: CPU %.2f s, Wall: %.2f s" % (cputime(c), walltime(w),)
            else:
                End = EndomorphismData(gtry, prec = digits)
            geo = End.geometric_representations()
            K = End.field_of_definition()
            g = gtry;
            break;
        except TypeError:
            pass

    if label == '810.a.196830.1':
	power += 2


    buffer =  "#working with the model y^2 = g(x)\n";
    buffer += "curve_dict['%s'] = %s;\n" % ('g', g.list());
    curve_dict['g'] = g.list();

    out += buffer;
    if verbose:
        print buffer

    # initialize numerical methods
    
    iaj = InvertAJglobal(g, prec, internalverbose);
    #CCap = iaj.C;

    buffer += "curve_dict['%s'] = %s;\n" % ('CCap', ' ComplexField(%s)' % curve_dict['prec']);
    #buffer += "CCap = curve_dict['CCap'];\n\n"

    # generator is r
    K = End.field_of_definition()
    if K.degree() == 1:
        K = QQ


    buffer = "# Field of definition of alpha\n"
    buffer += "curve_dict['%s'] = %s;\n" % ('K', sage_str_numberfield(K, 'x','r'));
    buffer += "K = %s\n" % ( "curve_dict['K']" )
    if K is not QQ:
        buffer += "r_approx = %s\n" % (K.gen().complex_embedding(),)
    out += buffer;
    out += "r = K.gen();\n"

    if verbose:
        print buffer




    #get the matrices over K
    geo = End.geometric_representations()

    #convert to sage
    # we also must take the transpose
    alphas = [  convert_magma_matrix_to_sage(y, K).transpose() for y in geo[0] ]
    alphas_geo = [ y.sage() for y in geo[2] ];
    
    d = len(alphas);
    out += "curve_dict['alphas_K'] = [None] * %d\n\n" % d
    out += "curve_dict['alphas_geo'] = [None] * %d\n\n" %d
    for i, _ in enumerate(alphas):
        alpha = alphas[i];
        alpha_geo = alphas_geo[i];
        buffer = "curve_dict['alphas_K'][%d] = Matrix(K, %s);\n" % (i, alpha.rows(),)
        buffer += "curve_dict['alphas_geo'][%d] = Matrix(%s);\n\n" % (i, alpha_geo.rows(),)
        out += buffer
        if verbose:
            print buffer
    curve_dict['alphas_K'] = alphas
    curve_dict['alphas_geo'] = alphas_geo

    out += "# where we stored all the data, for each alpha\n";
    out += "curve_dict['data'] = [{} for _ in range(%d) ] ;" % d
    curve_dict['data']  = [{} for _ in range(d) ];

    for i in range(d):
        if alphas[i] == Matrix([[1,0],[0,1]]):
            curve_dict['data'][i] = None;
        else:
            if verbose:
                print "Computing alpha(P + P)"
                print "where alpha = Matrix(K, %s)\n" %  (alphas[i].rows(),)

            # Pick a rational point, or a point over a quadratic extension with small discriminant
            for j, P0 in enumerate(find_rational_point(g)):
                if label in ['540800.a.540800.1', '529.a.529.1', '1521.a.41067.1','12500.a.12500.1', '18225.c.164025.1' ] and j < 2:
                    pass;
                elif label in ['810.a.196830.1'] and i == 2 and j == 0:
                    pass;
                else:
                    if verbose:
                        print "j = %s" % j
                    sys.stdout.flush()
                    sys.stderr.flush()
                    try:

                        output_alpha = compute_alpha_point(g, iaj, alphas[i], P0, digits, power, verbose = verbose, aggressive = True, append = str(i));

                        verified, trace_and_norm = add_trace_and_norm_ladic(g, output_alpha, alphas_geo[i], verbose = verbose);
                        break;
                    except ZeroDivisionError:
                        pass;
                    except OverflowError:
                        if verbose:
                            print "the mesh for numerical integral is too big"
                        pass;
                    if verbose:
                        print "trying with a new point\n"

            if verbose:
                print "\n\n\n"
            
            output_alpha['trace_and_norm'] = trace_and_norm;
            output_alpha['verified'] = verified
            output_alpha['alphas_geo'] = alphas_geo[i];
            internal_out = "\n"; #local_dict = {}\n\n";
            #
            
            
            # the fields
            internal_out += "# L the field where P and the algx_poly are defined\n"
            internal_out += "curve_dict['data'][%d]['L'] = %s\n" % (i, sage_str_numberfield(output_alpha['L'],'x','b'+str(i)),)
            if output_alpha['L'] is not QQ:
                internal_out += "b%d = curve_dict['data'][%d]['L'].gen()\n" % (i, i);
                internal_out += "curve_dict['data'][%d]['Lgen_approx'] = %s\n" % (i, output_alpha['L'].gen().complex_embedding());
            internal_out += "\n";

            internal_out += "curve_dict['data'][%d]['alpha'] = Matrix(curve_dict['data'][%d]['L'], %s) \n\n" % (i, i, [[elt for elt in row] for row in output_alpha['alpha'].rows()]);
            internal_out += "curve_dict['data'][%d]['alpha_geo'] = curve_dict['alphas_geo'][%d]\n\n" % (i, i)



            internal_out += "curve_dict['data'][%d]['P'] = vector(curve_dict['data'][%d]['L'], %s)\n\n" % (i, i, output_alpha['P']);

            
            
            internal_out += "# alpha(P + P) - \inf = R0 + R1 - \inf\n";
            
            internal_out +="\n\n\n";
            for key in ['algx_poly','R','x_poly', 'trace_and_norm', 'verified']:
                internal_out += "curve_dict['data'][%d]['%s'] = %s;\n" % (i, key, output_alpha[key]);
            
            # this just clutters the file, 
            # internal_out += "curve_dict['data'][%d]['ajP'] = vector(%s)\n\n" % (i, output_alpha['ajP']);
            # internal_out += "curve_dict['data'][%d]['alpha2P0'] = alphas_ap[%d] * 2 * local_dict['ajP']\n\n" % (i, i)

            internal_out += "\n\n\n"
            curve_dict['data'] [i] = output_alpha;
            out += internal_out;
    verified = all([elt['verified'] for elt in curve_dict['data'] if elt is not None])
    out += "curve_dict['verified'] = %s\n\n" % verified
    return curve_dict, out, verified;









# old stuff
#def add_newton_lift_ap(curve_dict, method = "ladic", verbose = True):
#    assert method == "ladic" or method == "ap"
#    
#    print "label = %s" % label
#    out = ""
#    for i, D in enumerate(data):
#        if D is not None:
#            if method == "ladic":
#                verified, trace_and_norm = add_trace_and_norm_ladic(D, verbose = verbose)
#            elif method == "ap":
#                verified, trace_and_norm = add_trace_and_norm_ap(D, verbose = verbose)
#
#            D['trace_numerator'], D['trace_denominator'], D['norm_numerator'], D['norm_denominator'] = trace_and_norm
#    return data
#
#def append_newton_lift(filename, verbose = True):
#    append = add_newton_lift(filename, verbose)
#    newfilename = filename[:-5] + "_wlift.sage"
#    f = open(newfilename,'w');
#    fold = open(filename,'r');
#    f.write(fold.read());
#    fold.close();
#    f.write(append);
#    f.close();
#    print "Finished writing into: %s" % (newfilename,)
#
#
#def append_newton_lift_standalone(filename, verbose = True):
#    import sys
#    output_filename = filename[:-5] + "_wlift.sage"
#    stdout_filename = filename[:-5] + "_wlift.stdout"
#    stderr_filename = filename[:-5] + "_wlift.stderr"
#
#    sys.stdout = open(stdout_filename, 'w', int(0))
#    sys.stderr = open(stderr_filename, 'w')
#
#
#    time append = add_newton_lift(filename, verbose)
#
#    f = open(output_filename,'w');
#    fold = open(filename,'r');
#    f.write(fold.read());
#    fold.close();
#    f.write(append);
#    f.close();
#    print "Finished writing into: %s" % (output_filename,)
#
#    sys.stdout.flush()
#    os.fsync(sys.stdout.fileno())
#    sys.stderr.flush()
#    os.fsync(sys.stderr.fileno())
#    sys.stdout = sys.__stdout__
#    sys.stderr = sys.__stderr__
#    print "Done: label = %s" % (label)
#    print os.popen("tail -n 2 %s | head -n 1" % stdout_filename).read()
#    return 0;
#









