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

inputfile = 'gce_genus3_hyperelliptic_endos.txt'
outputfile = 'special_curves_hyp.txt'

# Index of the criterion:
index = 2
# A priori boring output:
boring = ["[[-1,1],[[[-1,1],[[['I',[-1,1],1,1]],[1,-1],['RR'],'undef']]]]"]

interesting = 0
with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            linestrip = line.rstrip()
            linesplit = linestrip.split(':')
            if not linesplit[index] in boring:
                outputstream.write(line)
