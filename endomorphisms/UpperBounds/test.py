#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.

import os.path
from sage.all import PolynomialRing, ZZ, QQ
from endomorphisms.UpperBounds.utils import flatten_list, get_frob_list_HyperellipticCurve
from endomorphisms import endomorphisms_upper_bound


def fetch_gce_list(filename):
    data = open(filename, "r");
    output = {};
    QQx = PolynomialRing(QQ,'x')
    for line in data.readlines():
        linesplit =  line.rstrip('\n').split(":")
        if len(linesplit) == 3:
            int_id_str, fh_str, endo_str = linesplit
            int_id = int(int_id_str)
            fh = map(QQx, fh_str[1:-1].split(','))
            endo = eval(endo_str)
            if endo[0] == [-1,1]: # genus 3 format
                endo = endo[1];
            rr_alg = sorted(endo[0][1][2])
            output[int_id] = fh + [rr_alg]
    return output



database_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], "../../database/");

genus2_hyperelliptic = fetch_gce_list(os.path.join(database_path, "gce_genus2_hyperelliptic_endos.txt"));
# genus2_hyperelliptic = {169: [x**5 + x**4, x**3 + x + 1, ['M_2(RR)']],
# 249: [x**2 + x, x**3 + 1, ['RR']],
# 277: [x**5 - 9*x**4 + 14*x**3 - 19*x**2 + 11*x - 6, 1, ['RR']],
# 294: [x**4 + x**2, x**3 + 1, ['RR', 'RR']],
# 295: [x**5 - 40*x**3 + 22*x**2 + 389*x - 608, x**2 + x + 1, ['RR']],
# 349: [-x**3 - x**2, x**3 + x**2 + x + 1, ['RR']],
# 353: [x**2, x**3 + x + 1, ['RR']],
# 389: [x**5 + 2*x**4 + 2*x**3 + x**2, x + 1, ['RR']],
# 394: [2*x**5 + x**4 - 12*x**3 + 17*x - 9, x**3 + x, ['RR']]}

genus3_hyperelliptic = fetch_gce_list(os.path.join(database_path, "gce_genus3_hyperelliptic_endos.txt"));
# genus3_hyperelliptic = {3993: [x**7 + x**6 + x**5 + x**3 + x**2 + x, x**4 + x**2 + 1, ['M_2(RR)', 'RR']],
# 5911: [x**3 + x**2 + x, x**4 + x**3 + x**2 + 1, ['RR']],
# 7744: [x**5 - x**4 + x**3, x**4 + 1, ['RR', 'RR']],
# 7967: [x**7 + x**6 - x**4 - x**2 - x, x**4 + x**2 + x + 1, ['RR']],
# 8233: [x**7 - 8*x**5 - 4*x**4 + 18*x**3 - 3*x**2 - 16*x + 8,
#  x**4 + x**3 + x**2 + 1,
#  ['RR']],
# 8907: [x**3 + x**2, x**4 + x**3 + x + 1, ['RR']],
# 9449: [x**7 - 2*x**6 + 5*x**5 - 2*x**4 + x**3 + 4*x**2 + 2*x, x**2 + 1, ['RR']],
# 10487: [-x**6 + x**2 - x, x**4 + 1, ['RR']],
# 10704: [-x**7 - x**5 + x**2, x**4 + x**3 + x + 1, ['RR']],
# 11113: [-x**5 + x**3 + x**2, x**4 + x + 1, ['RR']],
# 11132: [-x**7 + 2*x**6 - 2*x**5 + 3*x**4 + x**3 + 2*x**2 + x,
#  x**4 + x**3 + x + 1,
#  ['RR', 'RR']],
# 11765: [-x**2 - x, x**4 + x**3 + x + 1, ['RR']]}

ZZx = PolynomialRing(ZZ, 'x');
x = ZZx.gen();
# the key has no special meaning
genus3_special = {0: [21*x**7 + 37506*x**5 + 933261*x**3 + 5841759*x, 0, ['CC', 'CC', 'CC']],
 1: [x**7 + 6*x**5 + 9*x**3 + x, 0, ['CC', 'CC', 'CC']],
 2: [16*x**7 + 357*x**5 - 819*x**3 + 448*x, 0, ['CC', 'CC', 'CC']],
 3: [-4*x**8 + 105*x**6 - 945*x**4 + 2100*x**2 - 5895*x + 420,
  x**4,
  ['CC', 'CC', 'CC']],
 4: [x**7 - 14*x**6 + 210*x**5 - 658*x**4 + 245*x**3 + 588*x**2 + 637*x - 686,
  0,
  ['RR', 'RR', 'RR']]};

genus2_special =  fetch_gce_list(os.path.join(database_path, "special_curves_hyp.txt"));


for i, examples in enumerate([genus2_hyperelliptic, genus3_hyperelliptic, genus3_special, genus2_special]):
    print("Batch %d" % i);
    print("Testing %d curves" % len(examples))
    for key, triple in examples.iteritems():
        f, h, RRendo = triple;
        frob_list = get_frob_list_HyperellipticCurve(f,h);
        endo_upper_bound = endomorphisms_upper_bound(frob_list);
        if endo_upper_bound[0] == False:
            print("%s: Failed to compute an upper bound!" % triple)
            print(endo_upper_bound[1])
            print("eta = %s, t = %s" % (endo_upper_bound[2], endo_upper_bound[3]))
            print("")
        else:
            RRcomputed = sorted(flatten_list([ elt[3] for elt in endo_upper_bound[4] ]))
            assert RRendo == RRcomputed, "Failed to match upper bound for (%s, %s), the computed upper bound was %s\n%s" % (key, triple, RRcomputed,endo_upper_bound)
            print( "%s : PASS" % (triple,))
    print("Batch %d: PASS! \n\n" % i);




