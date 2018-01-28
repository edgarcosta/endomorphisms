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

from sage.all import magma, QQ, PolynomialRing, NumberField
from Dictionaries import index_dictionary

def pretty_print_over_field_description(struct, g):
    return statements_all(struct, g)

def pretty_print_lattice_description(lat, g):
    _index_dict_ = index_dictionary()
    statements = [ ]
    base = lat[_index_dict_['base'] - 1]
    entries = lat[_index_dict_['entries'] - 1]
    for entry in entries:
        statement = ''
        field = entry[_index_dict_['field'] - 1]
        statement += "Over subfield %s:\n" % pretty_print_field(field, base)
        structure = entry[_index_dict_['structure'] - 1]
        statement += statements_all(structure, g)
        statements.append(statement)
    return '\n\n'.join(statements)

def pretty_print_polynomial_list(pol, base, str_var = 'x', str_gen = 'r'):
    if pol[0] in QQ:
        return str(PolynomialRing(QQ, str_var)(pol))
    R = PolynomialRing(QQ, 't')
    F = NumberField(R(base), str_gen)
    return str(PolynomialRing(F, str_var)(pol))

def pretty_print_ring(field, index):
    if field[0] in QQ:
        if len(field) == 2:
            return 'ZZ'
        elif len(field) == 3:
            c, b, a = field
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
    # Can fix base in the upcoming way because this function is only used
    # relative to QQ
    if index == 1:
        return 'Int (%s)' % pretty_print_field(field, [-1, 1])
    else:
        return 'Sub (%s, %s)' % (pretty_print_field(field, [-1, 1]), index)

def pretty_print_field(field, base, str_field = 'F'):
    if field[0] in QQ:
        base = [-1, 1]
        str_field = 'QQ'
        if len(field) == 3:
            c, b, a = field
            D = b**2 - 4*a*c
            return '%s (sqrt(%s))' % (str_field, D.squarefree_part())
    if len(field) == 2:
        return '%s' % str_field
    return '%s [x] / (%s)' % (str_field, pretty_print_polynomial_list(field, base))

def statements_all(desc, g):
    statements= [ ]
    statements.append(statement_endomorphisms_QQ(desc, g) + ' ' + statement_cm(desc, g))
    statements.append(statement_endomorphisms_ZZ(desc, g) + ' ' + statement_eichler(desc, g))
    statements.append(statement_endomorphisms_RR(desc, g))
    statements.append(statement_sato_tate_group(desc, g))
    statements.append(statement_gl2(desc, g) + '; ' + statement_simple(desc, g))
    return '\n'.join(statements)

def statement_endomorphisms_QQ(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    factors_QQ = desc[_index_dict_['factors_QQ']]
    statements = [ statement_factor_QQ(factor_QQ) for factor_QQ in factors_QQ ]
    statement =  ' x '.join(statements)
    return "End (J_%s) ox QQ: " % str_field + statement

def statement_endomorphisms_ZZ(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    factors_QQ = desc[_index_dict_['factors_QQ']]
    desc_ZZ = desc[_index_dict_['desc_ZZ']]
    if desc_ZZ[_index_dict_['index']] == 1:
        statements = [ statement_factor_ZZ_maximal(factor_QQ, desc_ZZ) for factor_QQ in factors_QQ ]
        statement = ' x '.join(statements)
    else:
        statement = statement_factors_ZZ_index(factors_QQ, desc_ZZ)
    return "End (J_%s):       " % str_field + statement

def statement_endomorphisms_RR(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    return "End (J_%s) ox RR: %s" % (str_field, ' x '.join(desc[_index_dict_['desc_RR']]))

def statement_factor_QQ(factor_QQ):
    _index_dict_ = index_dictionary()
    # TODO: Assumes g <= 3
    dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
    disc = factor_QQ[_index_dict_['disc']]
    albert_type = factor_QQ[_index_dict_['albert_type']]
    str_center = pretty_print_field(factor_QQ[_index_dict_['center']], [-1, 1])
    str_dim_sqrt = str(dim_sqrt)
    str_disc = str(disc)

    if albert_type == 'I':
        if dim_sqrt == 1:
            statement = str_center
        else:
            statement = "M_%s (%s)" % (str_dim_sqrt, str_center)

    elif albert_type == 'II':
        if dim_sqrt == 2:
            statement = "IndefQuat (%s, %s)"  % (str_center, str_disc)
        else:
            statement = "IndefAlg_%s (%s, %s)"  % (str_dim_sqrt, str_center, str_disc)

    elif albert_type == 'III':
        if dim_sqrt == 2:
            statement = "DefQuat (%s, %s)"  % (str_center, str_disc)
        else:
            statement = "DefAlg_%s (%s, %s)"  % (str_dim_sqrt, str_center, str_disc)

    elif albert_type == 'IV':
        if dim_sqrt == 1:
            statement = str_center
        elif disc == 1:
            statement = "M_%s (%s)" % (str_dim_sqrt, str_center)
        elif dim_sqrt == 2:
            statement = "Quat (%s, %s)"  % (str_center, str_disc)
        else:
            statement = "Alg_%s (%s, %s)"  % (str_dim_sqrt, str_center, str_disc)

    return statement

def statement_factor_ZZ_maximal(factor_QQ, desc_ZZ, str_field = 'K'):
    _index_dict_ = index_dictionary()
    # TODO: Assumes g <= 3
    albert_type = factor_QQ[_index_dict_['albert_type']]
    field_pretty = pretty_print_field(factor_QQ[_index_dict_['center']], [-1, 1])
    ring_pretty = pretty_print_ring(factor_QQ[_index_dict_['center']], 1)
    disc = factor_QQ[_index_dict_['disc']]
    dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
    if albert_type == 'I' or albert_type == 'IV':
        if dim_sqrt == 1:
            return ring_pretty
        # TODO: Next line can be done in greater generality over PIDs, but to
        # found out whether we are in such a situation costs time, so omitted
        # for now.
        elif (disc == 1) and field_pretty == 'QQ':
            return 'M_%s (%s)' % (dim_sqrt, ring_pretty)
    return 'Max (%s)' % statement_factor_QQ(factor_QQ)

def statement_factors_ZZ_index(factors_QQ, desc_ZZ, str_field = 'K'):
    _index_dict_ = index_dictionary()
    # TODO: Assumes g <= 3
    index = desc_ZZ[_index_dict_['index']]
    if len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        albert_type = factor_QQ[_index_dict_['albert_type']]
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        if albert_type == 'I' and dim_sqrt == 1:
            desc_field = factor_QQ[_index_dict_['center']]
            return pretty_print_ring(desc_field, index)
    return "Sub (End (J_%s) ox QQ, %s)" % (str_field, index)

def statement_sato_tate_group(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    sato_tate = desc[_index_dict_['sato_tate']]
    if sato_tate == "" or sato_tate == "undef":
        sato_tate = "not classified yet"
    return "Sato-Tate group: %s" % sato_tate

def statement_cm(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    factors_QQ = desc[_index_dict_['factors_QQ']]
    dimsum = 0
    for factor_QQ in factors_QQ:
        albert_type = factor_QQ[_index_dict_['albert_type']]
        dim_base = len(factor_QQ[_index_dict_['center']]) - 1
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        if albert_type == 'IV':
            dimsum +=  dim_base * dim_sqrt
    if (dimsum // 2) == g:
        return "(CM)"
    return ""

def statement_eichler(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
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

def statement_gl2(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    factors_QQ = desc[_index_dict_['factors_QQ']]
    dimsum = 0
    for factor_QQ in factors_QQ:
        dim_base = len(factor_QQ[_index_dict_['center']]) - 1
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        dimsum += dim_base * dim_sqrt
    if dimsum == g:
        return "of GL_2-type"
    return "not of GL_2-type"

def statement_simple(desc, g, str_field = 'K'):
    _index_dict_ = index_dictionary()
    factors_QQ = desc[_index_dict_['factors_QQ']]
    if len(factors_QQ) == 1:
        factor_QQ = factors_QQ[0]
        dim_sqrt = factor_QQ[_index_dict_['dim_sqrt']]
        disc = factor_QQ[_index_dict_['disc']]
        if dim_sqrt == 1 or disc != 1:
            return "simple"
    return "not simple"
