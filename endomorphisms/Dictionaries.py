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

from sage.all import magma

def index_dictionary():
    dikt = dict()

    # Magma indices for lattice
    dikt['base'] = 1
    dikt['entries'] = 2

    # Magma indices for lattice entries
    dikt['field'] = 1
    dikt['structure'] = 2

    # Magma indices for base and field key
    dikt['seq'] = 1
    dikt['magma'] = 2

    # Magma indices for OverField
    dikt['representation'] = 1
    dikt['algebra'] = 2
    dikt['description'] = 3

    # Magma indices for representation key
    dikt['tangent'] = 1
    dikt['homology'] = 2

    # Magma indices for algebra key
    dikt['alg_QQ'] = 1
    dikt['alg_ZZ'] = 2
    dikt['alg_RR'] = 3
    dikt['alg_ST'] = 4

    # Sage indices for description key (and subcases, see below)
    dikt['factors_QQ'] = 0
    dikt['desc_ZZ'] = 1
    dikt['desc_RR'] = 2
    dikt['sato_tate'] = 3

    # Sage indices for a factor_QQ
    dikt['albert_type'] = 0
    dikt['center'] = 1
    dikt['d'] = 2
    dikt['disc'] = 3
    dikt['m'] = 4

    # Sage indices for desc_ZZ key
    dikt['index'] = 0
    dikt['is_eichler'] = 1

    # Magma indices for decomposition.
    dikt['field'] = 1
    dikt['idem'] = 2
    dikt['factor'] = 3
    dikt['proj'] = 4

    # Sage indices for decomposition factor
    dikt['fac_type'] = 0
    dikt['fac_field'] = 1
    dikt['fac_coeffs'] = 2
    return dikt

def sagify_description(desc_list):
    return eval(repr(magma.SagifyDescription(desc_list)))

def dict_lattice(lattice):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['base'] = dict_base(lattice[_index_dict_['base']])
    dikt['entries'] = dict_entries(lattice[_index_dict_['entries']])
    return dikt

def desc_lattice(lattice):
    _index_dict_ = index_dictionary()
    desc = [ ]
    desc.append(desc_base(lattice[_index_dict_['base']]))
    desc.append(desc_entries(lattice[_index_dict_['entries']]))
    return desc

def dict_entries(entries):
    _index_dict_ = index_dictionary()
    dicts = [ ]
    for tup in entries:
        dikt = dict()
        dikt['field'] = dict_field(tup[_index_dict_['field']])
        dikt['structure'] = dict_structure(tup[_index_dict_['structure']])
        dicts.append(dikt)
    return dicts

def desc_sthash(sthash):
    return sagify_description(sthash)

def desc_entries(entries):
    _index_dict_ = index_dictionary()
    descs = [ ]
    for tup in entries:
        desc = [ ]
        desc.append(desc_field(tup[_index_dict_['field']]))
        desc.append(desc_structure(tup[_index_dict_['structure']]))
        descs.append(desc)
    return descs

def dict_structure(structure):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['representation'] = dict_rep(structure[_index_dict_['representation']])
    dikt['algebra'] = dict_alg(structure[_index_dict_['algebra']])
    desc = sagify_description(structure[_index_dict_['description']])
    dikt['description'] = dict_description(desc)
    return dikt

def desc_structure(structure):
    _index_dict_ = index_dictionary()
    return sagify_description(structure[_index_dict_['description']])

def dict_base(base):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['seq'] = base[_index_dict_['seq']]
    dikt['magma'] = base[_index_dict_['magma']]
    return dikt

def desc_base(base):
    _index_dict_ = index_dictionary()
    return sagify_description(base[_index_dict_['seq']])

def dict_field(field):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['seq'] = field[_index_dict_['seq']]
    dikt['magma'] = field[_index_dict_['magma']]
    return dikt

def desc_field(field):
    _index_dict_ = index_dictionary()
    return sagify_description(field[_index_dict_['seq']])

def dict_rep(rep):
    _index_dict_ = index_dictionary()
    return [ dict_gen(gen) for gen in rep ]

def dict_gen(gen):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['tangent'] = gen[_index_dict_['tangent']]
    dikt['homology'] = gen[_index_dict_['homology']]
    return dikt

def dict_alg(rep):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['alg_QQ'] = rep[_index_dict_['alg_QQ']]
    dikt['alg_ZZ'] = rep[_index_dict_['alg_ZZ']]
    dikt['alg_RR'] = rep[_index_dict_['alg_RR']]
    dikt['alg_ST'] = rep[_index_dict_['alg_ST']]
    return dikt

def dict_description(desc):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['factors_QQ'] = [ dict_factor_QQ(factor_QQ) for factor_QQ in desc[_index_dict_['factors_QQ']] ]
    dikt['desc_ZZ'] = dict_desc_ZZ(desc[_index_dict_['desc_ZZ']])
    dikt['desc_RR'] = dict_desc_RR(desc[_index_dict_['desc_RR']])
    dikt['sato_tate'] = desc[_index_dict_['sato_tate']]
    return dikt

def dict_factor_QQ(factor_QQ):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['albert_type'] = factor_QQ[_index_dict_['albert_type']]
    dikt['center'] = factor_QQ[_index_dict_['center']]
    dikt['d'] = factor_QQ[_index_dict_['d']]
    dikt['disc'] = factor_QQ[_index_dict_['disc']]
    dikt['m'] = factor_QQ[_index_dict_['m']]
    return dikt

def dict_desc_ZZ(desc_ZZ):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['index'] = desc_ZZ[_index_dict_['index']]
    dikt['is_eichler'] = desc_ZZ[_index_dict_['is_eichler']]
    return dikt

def dict_desc_RR(desc_RR):
    _index_dict_ = index_dictionary()
    dikt = dict()
    dikt['factors'] = desc_RR
    return dikt
