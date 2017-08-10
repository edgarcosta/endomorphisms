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
base_string = 'gce_genus3_nonhyperelliptic'
inputfile = base_string + '_endos.txt'

# Index of the representations:
index = 2

counter = 0
interesting = 0
with open(inputfile) as inputstream:
    for line in inputstream:
        counter += 1
        linestrip = line.rstrip()
        linesplit = linestrip.split(':')
        if linesplit[index] != "[[-1,1],[[[-1,1],[[['I',[-1,1],1,1]],[1,-1],['RR'],'undef']]]]":
            L = eval(linesplit[index])
            #if len(L[1][0][1][0]) == 1:
            if True:
                interesting += 1
                print linesplit[0]
                print linesplit[1]
                print pretty_print_lattice_description(eval(linesplit[index]), g = 3)
print interesting
