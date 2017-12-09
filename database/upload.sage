"""
 *  Upload to LMFDB
"""

# Password for LMFDB (not like I will write this here for real)
passwd = 'bla'

inputfile = 'gce_genus3_hyperelliptic_endos.txt'
# Index of label and endomorphism lattice:
label = 0
lattice = 3

geom = 0
base = -1
entries = _index_dict_['entries'] - 1
field = _index_dict_['field'] - 1
structure = _index_dict_['structure'] - 1
seq = _index_dict_['seq'] - 1
center = _index_dict_['center']
factors_QQ = _index_dict_['factors_QQ']
desc_ZZ = _index_dict_['desc_ZZ']
desc_RR = _index_dict_['desc_RR']
disc = _index_dict_['disc']
fac_type = _index_dict_['fac_type']
fac_field = _index_dict_['fac_field']
fac_coeffs = _index_dict_['fac_coeffs']

# Ambient needed in what follows:
R.<x> = PolynomialRing(QQ)
counter = 0
while counter_done < num_lines:
    import pymongo
    lmfdb = pymongo.MongoClient(host = 'localhost', port = int(37010))
    lmfdb.numberfields.authenticate('editor', passwd)
    lmfdb.elliptic_curves.authenticate('editor', passwd)
    lmfdb.genus2_curves.authenticate('editor', passwd)

    nfdb = lmfdb.numberfields.fields
    ecdb_QQ = lmfdb.elliptic_curves.curves
    ecdb_nf = lmfdb.elliptic_curves.nfcurves
    endodb = lmfdb.genus2_curves.test_upload

    with open(inputfile) as inputstream:
        counter = 0
        for line in inputstream:
            counter += 1
            print counter
            linestrip = line.rstrip()
            linesplit = linestrip.split(':')

            # Starting new dictionary for the curve:
            endodata = {'label': linesplit[label]}

            # Input endomorphism lattice and prettify
            L = eval(linesplit[lattice])
            # Forget about the base since we are over QQ:
            EDs = L[entries]
            N = len(EDs)
            for i in range(N):
                ED = EDs[i]
                # Use strings in subfield description to track it down:
                coeffs = str(ED[field][seq]).replace(' ', '').replace('[','').replace(']', '')
                nf = nfdb.find_one({'coeffs': coeffs})
                # Prepend label if it exists:
                if nf:
                    ED[field].append(nf['label'])
                else:
                    ED[field].append('')
                for factorQQ in ED[structure][factors_QQ]:
                    # Prepend label for center if it exists:
                    coeffs = str(factorQQ[center]).replace(' ', '').replace('[','').replace(']', '')
                    nf = nfdb.find_one({'coeffs': coeffs})
                    if nf:
                        factorQQ[center].append(nf['label'])
                    else:
                        factorQQ[center].append('')
            # Lattice is as before, with labels (perhaps empty) prepended
            endodata['lattice'] = EDs

            # Field of definition, labels at last place:
            endodata['fod_coeffs'] = EDs[geom][field][seq]
            endodata['fod_label'] = EDs[geom][field][-1]

            # Information over QQ and QQbar:
            endodata['factors_QQ_base'] = EDs[base][structure][factors_QQ]
            endodata['factors_QQ_geom'] = EDs[geom][structure][factors_QQ]
            endodata['desc_ZZ_base'] = EDs[base][structure][desc_ZZ]
            endodata['desc_ZZ_geom'] = EDs[geom][structure][desc_ZZ]
            endodata['desc_RR_base'] = EDs[base][structure][desc_RR]
            endodata['desc_RR_geom'] = EDs[geom][structure][desc_RR]

            # Simple if there is one factor of the algebra, not a quaternion
            # algebra:
            endodata['is_simple_base'] = (len(EDs[base][structure][factors_QQ]) == 1 \
                and EDs[base][structure][factorsQQ][0][disc] == 1)
            endodata['is_simple_geom'] = (len(EDs[geom][structure][factors_QQ]) == 1 \
                and EDs[geom][structure][factorsQQ][0][disc] == 1)

            # Splitting information:
            fac_descs = sage_eval(linesplit[decomp])
            fac_dicts = [ ]
            for fac_desc in fac_descs:

                # TODO: Only elliptic curves so far; modify this for the genus 2 case
                if fac_desc[fac_type] == 'ell':
                    fac_dict = dic()
                    # Convert field to list of integers:
                    fac_dict['fod_coeffs'] = [ int(c) for c in fac_desc[fac_field] ]
                    # Convert two lists of fractions to two lists of strings:
                    fac_dict['fac_coeffs'] = [ [ repr(c) for c in coeff ] for coeff in fac_desc[fac_coeffs] ]

                    # Try to find a label for the field:
                    coeffs = str(fac_desc[fac_coeffs]).replace(' ', '').replace('[','').replace(']', '')
                    nf = nfdb.find_one({'coeffs': coeffs})
                    field_label = ''
                    if nf:
                        field_label = nf['label']
                    fac_dict['fod_label'] = field_label

                    # If such a label is available, try to get the factors:
                    K.<r> = NumberField(fac_desc[fac_field])
                    E0 = EllipticCurve([ K(c) for c in fac_desc[fac_coeffs] ])

                    # Now look in the database for curves with the same
                    # j-invariant and conductor (norm) and check if they
                    # are isomorphic
                    # Over QQ:
                    if K.degree() == 1:
                        conductor = int(E0.conductor())
                        fac_dict['fac_condnorm'] = int(conductor)
                        jinv = repr(E0.j_invariant())
                        ecs = ecdb_QQ.find({'conductor': conductor, 'jinv': jinv})
                        fac_dict['fac_label'] = ''
                        for ec in ecs:
                            E = EllipticCurve([ QQ(a.encode('ascii')) for a in ec['ainvs'] ])
                            if E0.is_isomorphic(E):
                                fac_dict['fac_label'] = ec['lmfdb_label']

                    # Over a general field:
                    else:
                        conductor_norm = int(E0.conductor().norm())
                        fac_dict['condnorm'] = conductor_norm
                        jinv = [ repr(c) for c in E0.j_invariant().list() ]
                        ecs = ecdb_nf.find({'field_label': field_label, 'conductor_norm': conductor_norm, 'jinv': jinv})
                        for ec in ecs:
                            E = EllipticCurve([ K([ sage_eval(c) for c in a ]) for a in ec['ainvs'] ])
                            fac_dict['fac_label'] = ''
                            if E0.is_isomorphic(E):
                                fac_dict['fac_label'] = ec['label']

                fac_dicts.append(fac_dict)
            endodata['fac_dicts'] = fac_dicts

            # Finally we add the darn thing
            if endodb.find_one({'label': endodata['label']}):
                endodb.find_one_and_replace({'label': endodata['label']}, endodata)
            else:
                endodb.insert(endodata)

# Check at end:
endodb = lmfdb.genus2_curves.test_upload
print endodb.find_one({'label': '169.a.169.1'})
