# Inspects endomorphism representations in a list by pretty-printing a dummy
import os, shutil

inputfile = 'gce_genus3_nonhyperelliptic_special_endos.txt'

index = 2
sthashes = [ ]
with open(inputfile) as inputstream:
    for line in inputstream:
        linestrip = line.rstrip()
        linesplit = linestrip.split(':')
        L = eval(linesplit[index])
        sthashes.append(tuple([ tuple(L[2][0]), tuple(L[2][1]) ]))

import collections
counter=collections.Counter(sthashes)
print(counter)

index = 2
sthashes = [ ]
with open(inputfile) as inputstream:
    for line in inputstream:
        linestrip = line.rstrip()
        linesplit = linestrip.split(':')
        L = eval(linesplit[index])
        gp = tuple(L[2][0])
        comp = tuple(sorted(L[1][0][1][2]))
        sig = tuple(L[2][1])
        sthashes.append(tuple([ gp, comp, sig ]))

import collections
counter=collections.Counter(sthashes)
print(counter)

triples = counter.keys()
for pair in cartesian_product([triples, triples]):
    trip1 = pair[0]
    trip2 = pair[1]
    if trip1[0] == trip2[0] and trip1[1] == trip2[1] and trip1[2] != trip2[2]:
        print "-------"
        print trip1, counter[trip1]
        print trip2, counter[trip1]

index = 2
sthashes = [ ]
with open(inputfile) as inputstream:
    for line in inputstream:
        linestrip = line.rstrip()
        linesplit = linestrip.split(':')
        L = eval(linesplit[index])
        gp = tuple(L[2][0])
        sig = tuple(L[2][1])
        comp = tuple(sorted(L[1][0][1][2]))
        comps = tuple(sorted([ tuple(sorted(c[1][2])) for c in L[1] ]))
        sthashes.append(tuple([ gp, sig, comps ]))

import collections
counter=collections.Counter(sthashes)
print(counter)

triples = counter.keys()
for pair in cartesian_product([triples, triples]):
    trip1 = pair[0]
    trip2 = pair[1]
    if trip1[0] == trip2[0] and trip1[1] == trip2[1] and trip1[2] != trip2[2]:
        print "-------"
        print trip1, counter[trip1]
        print trip2, counter[trip1]

