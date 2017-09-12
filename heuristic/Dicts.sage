"""
 *  Dictionary conversions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

_index_dict_ = dict()

# Magma indices for lattice
_index_dict_['base'] = 1
_index_dict_['entries'] = 2

# Magma indices for lattice entries
_index_dict_['field'] = 1
_index_dict_['structure'] = 2

# Magma indices for base and field key
_index_dict_['seq'] = 1
_index_dict_['magma'] = 2
# Class number not used because of performance:
#_index_dict_['h'] = 3

# Magma indices for OverField
_index_dict_['representation'] = 1
_index_dict_['algebra'] = 2
_index_dict_['description'] = 3

# Magma indices for representation key
_index_dict_['tangent'] = 1
_index_dict_['homology'] = 2
_index_dict_['approx'] = 3

# Magma indices for algebra key
_index_dict_['alg_QQ'] = 1
_index_dict_['alg_ZZ'] = 2
_index_dict_['alg_RR'] = 3
_index_dict_['alg_ST'] = 4

# Sage indices for description key (and subcases, see below)
_index_dict_['factors_QQ'] = 0
_index_dict_['desc_ZZ'] = 1
_index_dict_['desc_RR'] = 2
_index_dict_['sato_tate'] = 3

# Sage indices for a factor_QQ
_index_dict_['albert_type'] = 0
_index_dict_['center'] = 1
_index_dict_['dim_sqrt'] = 2
_index_dict_['disc'] = 3

# Sage indices for desc_ZZ key
_index_dict_['index'] = 0
_index_dict_['is_eichler'] = 1

# Magma indices for decomposition.
_index_dict_['field'] = 1
_index_dict_['idem'] = 2
_index_dict_['factor'] = 3
_index_dict_['proj'] = 4

# Sage indices for decomposition factor
_index_dict_['fac_type'] = 0
_index_dict_['fac_field'] = 1
_index_dict_['fac_coeffs'] = 2

def sagify_description(desc_list):
    return eval(repr(magma.SagifyDescription(desc_list)))

def dict_lattice(lattice):
    dikt = dict()
    dikt['base'] = dict_base(lattice[_index_dict_['base']])
    dikt['entries'] = dict_entries(lattice[_index_dict_['entries']])
    return dikt

def desc_lattice(lattice):
    desc = [ ]
    desc.append(desc_base(lattice[_index_dict_['base']]))
    desc.append(desc_entries(lattice[_index_dict_['entries']]))
    return desc

def dict_entries(entries):
    dicts = [ ]
    for tup in entries:
        dikt = dict()
        dikt['field'] = dict_field(tup[_index_dict_['field']])
        dikt['structure'] = dict_structure(tup[_index_dict_['structure']])
        dicts.append(dikt)
    return dicts

def desc_entries(entries):
    descs = [ ]
    for tup in entries:
        desc = [ ]
        desc.append(desc_field(tup[_index_dict_['field']]))
        desc.append(desc_structure(tup[_index_dict_['structure']]))
        descs.append(desc)
    return descs

def dict_structure(structure):
    dikt = dict()
    dikt['representation'] = dict_rep(structure[_index_dict_['representation']])
    dikt['algebra'] = dict_alg(structure[_index_dict_['algebra']])
    desc = sagify_description(structure[_index_dict_['description']])
    dikt['description'] = dict_description(desc)
    return dikt

def desc_structure(structure):
    return sagify_description(structure[_index_dict_['description']])

def dict_base(base):
    dikt = dict()
    dikt['seq'] = base[_index_dict_['seq']]
    dikt['magma'] = base[_index_dict_['magma']]
    return dikt

def desc_base(base):
    return sagify_description(base[_index_dict_['seq']])

def dict_field(field):
    dikt = dict()
    dikt['seq'] = field[_index_dict_['seq']]
    dikt['magma'] = field[_index_dict_['magma']]
    return dikt

def desc_field(field):
    return sagify_description(field[_index_dict_['seq']])

def dict_rep(rep):
    return [ dict_gen(gen) for gen in rep ]

def dict_gen(gen):
    dikt = dict()
    dikt['tangent'] = gen[_index_dict_['tangent']]
    dikt['homology'] = gen[_index_dict_['homology']]
    dikt['approx'] = gen[_index_dict_['approx']]
    return dikt

def dict_alg(rep):
    dikt = dict()
    dikt['alg_QQ'] = rep[_index_dict_['alg_QQ']]
    dikt['alg_ZZ'] = rep[_index_dict_['alg_ZZ']]
    dikt['alg_RR'] = rep[_index_dict_['alg_RR']]
    dikt['alg_ST'] = rep[_index_dict_['alg_ST']]
    return dikt

def dict_description(desc):
    dikt = dict()
    dikt['factors_QQ'] = [ dict_factor_QQ(factor_QQ) for factor_QQ in desc[_index_dict_['factors_QQ']] ]
    dikt['desc_ZZ'] = dict_desc_ZZ(desc[_index_dict_['desc_ZZ']])
    dikt['desc_RR'] = dict_desc_RR(desc[_index_dict_['desc_RR']])
    dikt['sato_tate'] = desc[_index_dict_['sato_tate']]
    return dikt

def dict_factor_QQ(factor_QQ):
    dikt = dict()
    dikt['albert_type'] = factor_QQ[_index_dict_['albert_type']]
    dikt['center'] = factor_QQ[_index_dict_['center']]
    dikt['dim_sqrt'] = factor_QQ[_index_dict_['dim_sqrt']]
    dikt['disc'] = factor_QQ[_index_dict_['disc']]
    return dikt

def dict_desc_ZZ(desc_ZZ):
    dikt = dict()
    dikt['index'] = desc_ZZ[_index_dict_['index']]
    dikt['is_eichler'] = desc_ZZ[_index_dict_['is_eichler']]
    return dikt

def dict_desc_RR(desc_RR):
    dikt = dict()
    dikt['factors'] = desc_RR
    return dikt
