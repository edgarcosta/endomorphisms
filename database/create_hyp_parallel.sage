# Adds endomorphism data to a file of colon-separated lines
from endomorphisms import *

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Length of the lines in the input file:
line_length = 3
# Specify indices of defining polynomials and Sato-Tate group here;
# making the latter negative ignores the corresponding check.
fh_index = 2
st_index = -1
# Precision:
prec = 300

import os, shutil

# Specify input and output:
base_string = 'gce_genus3_hyperelliptic_possibly_special'
inputfile = base_string + '.txt'
outputfile = base_string + '_endos.txt'

# Ambient ring:
R.<x> = PolynomialRing(QQ)

@parallel(ncpus=4)
def step(line):
    global fh_index
    global R
    global prec
    global outputstream
    try:
        linestrip = line.rstrip()
        print linestrip
        linesplit = linestrip.split(':')
        pol_list = eval(linesplit[fh_index].replace('^', '**'))
        f = R(pol_list[0])
        h = R(pol_list[1])
        X = HyperellipticCurve(f, h)
        Endo = EndomorphismData(X, prec = prec, molin_neurohr = True)
        lat_str = Endo.lattice()._desc_
        sth_str = Endo.lattice()._sthash_
        line_new1 = repr(lat_str).replace('\n', '').replace(' ', '')
        line_new2 = repr(sth_str).replace('\n', '').replace(' ', '')
        return line_new1 + ':' + line_new2
    except:
        print linestrip
        print "Error"
        return "Error"

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        inputlist = [ line for line in inputstream ]
        outputlist = [ ]
        for tup in list(step(inputlist)):
            outputlist.append(tup[0][0][0].rstrip() + ':' + tup[1] + '\n')
        outputlist.sort(key = lambda line: int(line.split(':')[0]))
        for line in outputlist:
            outputstream.write(line)
