# Inspects endomorphism representations in a list by pretty-printing a dummy
import os, shutil

inputfile = 'gce_genus3_nonhyperelliptic_endos.txt'

# Index of the representations:
index = 2
# A priori boring output:
boring = ["[[-1,1],[[[-1,1],[[['I',[-1,1],1,1]],[1,-1],['RR'],'undef']]]]"]

geom = 0
base = -1
entries = _index_dict_['entries'] - 1
structure = _index_dict_['structure'] - 1
factors_QQ = _index_dict_['factors_QQ']
desc_ZZ = _index_dict_['desc_ZZ']
desc_ZZ_index = _index_dict_['index']
disc = _index_dict_['disc']

interesting = 0
with open(inputfile) as inputstream:
    for line in inputstream:
        linestrip = line.rstrip()
        linesplit = linestrip.split(':')
        if not linesplit[index] in boring:
            L = eval(linesplit[index])
            EDs = L[entries]
            # All interesting cases:
            if True:
            # Geometrically simple:
            #if len(EDs[geom][structure][factors_QQ]) == 1 and EDs[geom][structure][factorsQQ][0][disc] == 1:
            # Three distinct geometric factors:
            #if len(EDs[geom][structure][factors_QQ]) == 3:
            # Two distinct geometric factors:
            #if len(EDs[geom][structure][factors_QQ]) == 2:
            # Two distinct geometric factors and index divisibility:
            #if (len(EDs[geom][structure][factors_QQ]) == 2) and (EDs[geom][structure][desc_ZZ][desc_ZZ_index] % 2 == 0):
                interesting += 1
                print ""
                print linesplit[0:index]
                print L
                print pretty_print_lattice_description(L, g = 3)
print ""
print "Total number of entries:"
print interesting
