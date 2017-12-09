# Inspects endomorphism representations in a list by pretty-printing a dummy
import os, shutil

inputfile = 'gce_genus3_nonhyperelliptic_endos.txt'
outputfile = 'special_curves_nonhyp.txt'

# Index of the criterion:
index = 1
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
