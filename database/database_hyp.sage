"""
 *  Generates endomorphism data from a colon-separated list.
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

# Adds endomorphism data to a file of colon-separated lines

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Length of the lines in the input file:
line_length = 2
# Specify indices of defining polynomials and Sato-Tate group here;
# making the latter negative ignores the corresponding check.
fh_index = 1
st_index = -1
# Precision:
prec = 300

import os, shutil

# Specify input and output:
base_string = 'gce_genus3_hyperelliptic'
inputfile = base_string + '.txt'
intermediatefile = base_string + '_temp.txt'
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
        linesplit = linestrip.split(':')
        pol_list = eval(linesplit[fh_index].replace('^', '**'))
        f = R(pol_list[0])
        h = R(pol_list[1])
        X = HyperellipticCurve(f, h)
        # Regular version:
        Endo = EndomorphismData(X, prec = prec, have_oldenburg = True)
        Lat_str = Endo.lattice()._desc_
        return repr(Lat_str).replace('\n', '').replace(' ', '')
    except:
        return "Error"

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        inputlist = [ line for line in inputstream ]
        for tup in list(step(inputlist)):
            outputstream.write(tup[0][0][0].rstrip() + ':' + tup[1] + '\n')
