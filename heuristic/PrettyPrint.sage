"""
 *  Pretty print functions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

def pretty_print_over_field_description(desc, genus, str_field):
    return statements_all(desc, genus, str_field)

def pretty_print_lattice_description(desc_lat, genus, str_field, str_var):
    statements = [ ]
    for desc in desc_lat:
        statement = ''
        desc_field = desc[_index_dict_['field'] - 1]
        statement += "Over subfield %s:\n" % pretty_print_field(desc_field)
        desc_structure = desc[_index_dict_['structure'] - 1]
        statement += statements_all(desc_structure, genus, str_field)
        statements.append(statement)
    return '\n\n'.join(statements)

def pretty_print_polynomial_list(s, str_var):
    return str(PolynomialRing(QQ, str_var)(s))

def pretty_print_ring(desc_field, index):
    # TODO: General base field, and for example cyclotomic extensions
    if len(desc_field) == 2:
        return 'ZZ'
    elif len(desc_field) == 3:
        c ,b, a = desc_field
        disc = b**2 - 4*a*c
        if index == 1:
            if disc % 4 == 0:
                return 'ZZ [sqrt(%s)]'% str(disc//4)
            return 'ZZ [(1 + sqrt(%s))/2]' % str(disc)
        if disc % 4 == 0:
            return 'ZZ [%s sqrt(%s)]' % (str(index), str(disc//4))
        if index % 2 == 0:
            return 'ZZ [%s sqrt(%s)]' % (str(index//2), str(disc))
        return 'ZZ [(1 + sqrt(%s))/2]' % (str(index), str(disc))
    else:
        if index == 1:
            return 'Int (%s)' % pretty_print_field(desc_field)
        else:
            return 'Sub (%s, %s)' % (pretty_print_field(desc_field), index)

def pretty_print_field(desc_field):
    # Only makes a difference for at most quadratic fields
    # TODO: General base field, and for example cyclotomic extensions
    if len(desc_field) == 2:
        return 'QQ'
    if len(desc_field) == 3:
        c, b, a = desc_field
        D = b**2 - 4*a*c
        return 'QQ (sqrt(%s))' % D.squarefree_part()
    else:
        return 'QQ [x] / (%s)' % pretty_print_polynomial_list(desc_field, 'x')

def statements_all(desc, genus, str_field):
    statements= [ ]
    statements.append(statement_endomorphisms_QQ(desc, genus, str_field) + ' ' + statement_cm(desc, genus, str_field))
    statements.append(statement_endomorphisms_ZZ(desc, genus, str_field) + ' ' + statement_eichler(desc, genus, str_field))
    statements.append(statement_endomorphisms_RR(desc, genus, str_field))
    statements.append(statement_sato_tate_group(desc, genus, str_field))
    statements.append(statement_gl2(desc, genus, str_field) + '; ' + statement_simple(desc, genus, str_field))
    return '\n'.join(statements)

def statement_endomorphisms_QQ(desc, genus, str_field):
    factor_dicts = desc[_index_dict_['factors_QQ']]
    statements = [ statement_factor_QQ(factor_dict) for factor_dict in factor_dicts ]
    statement =  ' x '.join(statements)
    return "End (J_%s) ox QQ: " % str_field + statement

def statement_endomorphisms_ZZ(desc, genus, str_field):
    factors_QQ = desc[_index_dict_['factors_QQ']]
    desc_ZZ = desc[_index_dict_['desc_ZZ']]
    if desc_ZZ[_index_dict_['index']] == 1:
        statements = [ statement_factor_ZZ_maximal(factor_QQ, desc_ZZ, str_field) for factor_QQ in factors_QQ ]
        statement = ' x '.join(statements)
    else:
        statement = statement_factors_ZZ_index(factors_QQ, desc_ZZ, str_field)
    return "End (J_%s):       " % str_field + statement

def statement_endomorphisms_RR(desc, genus, str_field):
    return "End (J_%s) ox RR: %s" % (str_field, ' x '.join(desc[_index_dict_['desc_RR']]))

def statement_factor_QQ(factor_QQ):
    # FIXME: Assumes g <= 3
    dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
    disc = factor_QQ[_index_dict_['disc']]
    albert_type = factor_QQ[_index_dict_['albert_type']]
    str_base = pretty_print_field(factor_QQ[_index_dict_['base_field']])
    str_dim_sqrt = str(dim_sqrt)
    str_disc = str(disc)

    if albert_type == 'I':
        if dim_sqrt == 1:
            statement = str_base
        else:
            statement = "M_%s (%s)" % (str_dim_sqrt, str_base)

    elif albert_type == 'II':
        if dim_sqrt == 2:
            statement = "IndefQuat (%s, %s)"  % (str_base, str_disc)
        else:
            statement = "IndefAlg_%s (%s, %s)"  % (str_dim_sqrt, str_base, str_disc)

    elif albert_type == 'III':
        if dim_sqrt == 2:
            statement = "DefQuat (%s, %s)"  % (str_base, str_disc)
        else:
            statement = "DefAlg_%s (%s, %s)"  % (str_dim_sqrt, str_base, str_disc)

    elif albert_type == 'IV':
        if dim_sqrt == 1:
            statement = str_base
        elif disc == 1:
            statement = "M_%s (%s)" % (str_dim_sqrt, str_base)
        elif dim_sqrt == 2:
            statement = "Quat (%s, %s)"  % (str_base, str_disc)
        else:
            statement = "Alg_%s (%s, %s)"  % (str_dim_sqrt, str_base, str_disc)

    return statement

def statement_factor_ZZ_maximal(factor_QQ, desc_ZZ, str_field):
    # FIXME: Assumes g <= 3
    albert_type = factor_QQ[_index_dict_['albert_type']]
    field_pretty = pretty_print_field(factor_QQ[_index_dict_['base_field']])
    ring_pretty = pretty_print_ring(factor_QQ[_index_dict_['base_field']], 1)
    disc = factor_QQ[_index_dict_['disc']]
    dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
    if albert_type == 'I' or albert_type == 'IV':
        if dim_sqrt == 1:
            return ring_pretty
        # TODO: Next line in greater generality over PIDs
        elif (disc == 1) and field_pretty == 'QQ':
            return 'M_%s (%s)' % (dim_sqrt, ring_pretty)
    return 'Max (%s)' % statement_factor_QQ(factor_QQ)

def statement_factors_ZZ_index(factors_QQ, desc_ZZ, str_field):
    # FIXME: Assumes g <= 3
    index = desc_ZZ[_index_dict_['index']]
    if len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        albert_type = factor_QQ[_index_dict_['albert_type']]
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        if albert_type == 'I' and dim_sqrt == 1:
            desc_field = factor_QQ[_index_dict_['base_field']]
            return pretty_print_ring(desc_field, index)
    return "Sub (End (J_%s) ox QQ, %s)" % (str_field, index)

def statement_sato_tate_group(desc, genus, str_field):
    sato_tate = desc[_index_dict_['sato_tate']]
    if sato_tate == "" or sato_tate == "undef":
        sato_tate = "not classified yet"
    return "Sato-Tate group: %s" % sato_tate

def statement_cm(desc, genus, str_field):
    factors_QQ = desc[_index_dict_['factors_QQ']]
    dimsum = 0
    for factor_QQ in factors_QQ:
        albert_type = factor_QQ[_index_dict_['albert_type']]
        dim_base = len(factor_QQ[_index_dict_['base_field']]) - 1
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        if albert_type == 'IV':
            dimsum +=  dim_base * dim_sqrt
    if (dimsum // 2) == genus:
        return "(CM)"
    return ""

def statement_eichler(desc, genus, str_field):
    factors_QQ = desc[_index_dict_['factors_QQ']]
    desc_ZZ = desc[_index_dict_['desc_ZZ']]
    if len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        albert_type = factor_QQ[_index_dict_['albert_type']]
        is_eichler = desc_ZZ[_index_dict_['is_eichler']]
        if is_eichler == 1:
            return "(Eichler)"
        elif is_eichler == 0:
            return "(non-Eichler)"
    return ""

def statement_gl2(desc, genus, str_field):
    factors_QQ = desc[_index_dict_['factors_QQ']]
    dimsum = 0
    for factor_QQ in factors_QQ:
        dim_base = len(factor_QQ[_index_dict_['base_field']]) - 1
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        dimsum += dim_base * dim_sqrt
    if dimsum == genus:
        return "of GL_2-type"
    return "not of GL_2-type"

def statement_simple(desc, genus, str_field):
    factors_QQ = desc[_index_dict_['factors_QQ']]
    if len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        disc = factor_QQ[_index_dict_['disc']]
        if dim_sqrt == 1 or disc != 1:
            return "simple"
    return "not simple"
