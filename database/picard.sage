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

# Inspects endomorphism representations in a list by pretty-printing a dummy
import os, shutil

# Specify input and output:
inputfile = 'gce_genus3_nonhyperelliptic_possibly_special.txt'
outputfile = 'picard.txt'

# Index of the curve:
index = 1

R.<x,y,z> = PolynomialRing(QQ, 3)
Rp.<x,y,z> = PolynomialRing(FiniteField(37), 3)
counter = 0
indices = [ 1, 2, 5, 6, 7, 8, 11, 12 ]
with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            counter += 1
            print counter
            linestrip = line.rstrip()
            linesplit = linestrip.split(':')
            pol_list = eval(linesplit[index].replace('^', '**'))
            Fp = magma(Rp(pol_list[0]))
            DOp = magma.DixmierOhnoInvariants(Fp)
            testp = all([ DOp[i] == 0 for i in indices ])
            if testp:
                F = magma(R(pol_list[0]))
                DO = magma.DixmierOhnoInvariants(F)
                test = all([ DO[i] == 0 for i in indices ])
                if test:
                    outputstream.write(line)
