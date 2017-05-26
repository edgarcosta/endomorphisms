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
# Precision (not below 200 please):
prec = 300

import os, shutil

# Specify input and output:
base_string = 'gce_genus3_hyperelliptic'
inputfile = base_string + '.txt'
intermediatefile = base_string + '_temp.txt'
outputfile = base_string + '_endos.txt'

# Ambient ring:
R.<x> = PolynomialRing(QQ)
# Default substitution (currently the identity only):
subst_list = [ x ]
# Bound on height of random substitutions:
B = 3
# The maximum number of runs:
maxrun = 1
# The counter for the current run and the line:
run = 0
counter = 0
done_list = [ ]

stop = False
exhaust = False
while not stop:
    run += 1
    counter = 0
    if not exhaust:
        subst = subst_list.pop()
    else:
        while True:
            num = R([QQ.random_element(B) for i in range (2)])
            den = R([QQ.random_element(B) for i in range (2)])
            if num != 0 and den != 0:
                break
            # We do not have to guard against the constant case, since that
            # bugs in Magma already and will get caught below.
        subst = num/den
    stop = True
    print "Run:", run
    print "Substitution:", subst
    with open(inputfile) as inputstream:
        with open(outputfile, 'w') as outputstream:
            for line in inputstream:
                counter += 1
                linestrip = line.rstrip()
                linesplit = linestrip.split(':')
                linestart = linestrip
                # We have to see if there is no new information on the line yet:
                if not counter in done_list:
                    print counter
                    pol_list = eval(linesplit[fh_index].replace('^', '**'))
                    f = R(pol_list[0])
                    h = R(pol_list[1])
                    den = subst.denominator()
                    f = R(den^6 * f(subst))
                    h = R(den^3 * h(subst))
                    X = HyperellipticCurve(f, h)
                    try:
                        # Regular version:
                        Endo = EndomorphismData(X, prec = prec, have_oldenburg = True)
                        Lat_str = Endo.lattice()._desc_
                        outputstream.write(linestart
                               + ':' + repr(Lat_str).replace('\n', '').replace(' ', '')
                               + '\n')
                        done_list.append(counter)
                        # Or check Sato-Tate and write if a match occurs:
                        #if st_index < 0 or Lat_str[-1][4] == linesplit[st_index]:
                        #    outputstream.write(linestart
                        #           + ':' + repr(Lat_str).replace('\n', '').replace(' ', '')
                        #           + '\n')
                        #    done_list.append(counter)
                        #else:
                        #    # In case of incorrect ST postpone until next time:
                        #    outputstream.write(line)
                        # Or make the database smaller:
                        #outputstream.write(linesplit[0] + ':' + linesplit[3] + '\n')
                        # Or check for lacking entries:
                        #if len(linesplit) != 3:
                        #    outputstream.write(line)
                        # Or check for large lattice
                        #if len(linesplit) == 3:
                        #    L = sage_eval(linesplit[2])
                        #    if len(L) != 1:
                        #        outputstream.write(line)
                    except:
                        # In case of an error postpone until next time:
                        print "Error"
                        outputstream.write(line)
                        stop = False
                else:
                    # Skip the line if it has been calculated already:
                    outputstream.write(line)
    inputfile = intermediatefile
    shutil.copyfile(outputfile, inputfile)
    if run >= maxrun:
        stop = True
os.remove(intermediatefile)
