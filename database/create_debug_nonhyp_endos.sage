# Adds endomorphism data to a file of colon-separated lines

# Defining polynomials have to be provided in pairs, defined by strings in x or
# by lists of integers. These polynomials (and the conjectural Sato-Tate group,
# if provided) need to be at a consistent index in the provided lines.

# Specify indices of defining polynomials
fh_index = 2
# Precision:
prec = 100

import os, shutil

# Specify input and output:
base_string = 'special_curves_nonhyp'
inputfile = base_string + '.txt'
outputfile = base_string + '_endos.txt'

# Ambient ring:
R.<x,y,z> = QQ[]

with open(inputfile) as inputstream:
    with open(outputfile, 'w') as outputstream:
        for line in inputstream:
            print line
            try:
                linestrip = line.rstrip()
                linesplit = linestrip.split(':')
                pol_list = eval(linesplit[fh_index].replace('^', '**'))
                f = R(pol_list[0])
                X = mPlaneCurve(f)
                Endo = EndomorphismData(X, prec = prec, have_oldenburg = True)
                test = Endo.verify()
                if not test:
                    print 'False'
                    outputstream.write(linestrip + ':' + 'False' + '\n')
            except:
                print 'Error'
                outputstream.write(linestrip + ':' + 'Error' + '\n')
