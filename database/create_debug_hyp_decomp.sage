# Adds endomorphism data to a file of colon-separated lines

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Specify indices of defining polynomials
fh_index = 1
# Precision:
prec = 300

import os, shutil

# Specify input and output:
base_string = 'special_curves_hyp'
inputfile = base_string + '.txt'
outputfile = base_string + '_decomp.txt'

# Ambient ring:
R.<x> = PolynomialRing(QQ)

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            print line
            try:
                linestrip = line.rstrip()
                linesplit = linestrip.split(':')
                pol_list = eval(linesplit[fh_index].replace('^', '**'))
                f = R(pol_list[0])
                h = R(pol_list[1])
                X = HyperellipticCurve(f, h)
                Endo = EndomorphismData(X, prec = prec, have_oldenburg = True)
                Dec = Endo.decomposition()
                facs = Dec.factors()
                idems = Dec.idempotents()
                test = Dec.verify()
                if not test:
                    print 'False'
                    outputstream.write(linestrip + ':' + 'False' + '\n')

            except:
                print 'Error'
                outputstream.write(linestrip + ':' + 'Error' + '\n')
