# Adds endomorphism data to a file of colon-separated lines
from endomorphisms import *

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Length of the lines in the input file:
line_length = 2
# Specify indices of defining polynomials and Sato-Tate group here;
# making the latter negative ignores the corresponding check.
f_index = 1
st_index = -1
# Precision:
prec = 230

import os, shutil

# Specify input and output:
base_string = 'gce_genus3_nonhyperelliptic_special'
inputfile = base_string + '_inter.txt'
outputfile = base_string + '_purged.txt'

# Ambient ring:
R.<x,y,z> = PolynomialRing(QQ, 3)

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            linestrip = line.rstrip()
            linesplit = linestrip.split(':')
            if len(linesplit) != line_length:
                outputstream.write(linesplit[0] + ':' + linesplit[1] + '\n')
