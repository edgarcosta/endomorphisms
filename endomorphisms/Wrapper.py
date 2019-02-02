"""
 *  Class wrappers for the algorithms
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

from Curves import *
from Dictionaries import *
from PrettyPrint import *
from Representations import *

from sage.all import magma

class EndomorphismData:
    def __init__(self, X, periods = ""):
        self.X = X
        self.g = magma.Genus(self.X)
        self.F = magma.BaseRing(self.X)
        if periods:
            self._P_ = periods
        else:
            self._P_ = self.period_matrix()
        self._calculate_geometric_representation_()

    def __repr__(self):
        return repr_endomorphism_data(self)

    def period_matrix(self):
        self._P_ = magma.PeriodMatrix(self.X)
        return self._P_

    def _calculate_geometric_representation_(self):
        self._geo_rep_list_ = magma.GeometricEndomorphismRepresentation(self._P_, self.F)
        self._geo_rep_dict_ = dict_rep(self._geo_rep_list_)

    def endomorphism_field(self):
        return magma.BaseRing(self._geo_rep_list_[1][1])

    def lattice(self):
        if not hasattr(self, "_lat_"):
            self._lat_ = Lattice(self)
        return self._lat_

    def decomposition(self):
        if not hasattr(self, "_lat_"):
            self._lat_ = self.lattice()
        return Decomposition(self)

#  TODO: Get this to work and uncomment it
#    def verify_algebra(self):
#        self._test_alg_ = True
#        return self._test_alg_
#
#    def verify_algebra_NS(self):
#        assert self.g == 2, "for now the upper bounds are only implemented for genus = 2";
#        # X: y^2 = g(x)
#        g = PolynomialRing(QQ, "x")(magma.Coefficients(magma.HyperellipticPolynomials(magma.SimplifiedModel(X))));
#        factorsRR_geom = self._desc_[self._index_dict_['desc_RR']];
#        algebra = self._list_[self._index_dict_['algebra']][self._index_dict_['alg_QQ']];
#        if factorsRR_geom ==  ['RR', 'RR'] and len(magma.DirectSumDecomposition(algebra)) == 1:
#            RM_coeff =  list(magma.Coefficients(magma.DefiningPolynomial(magma.BaseRing(magma.AlgebraOverCenter(algebra)))))
#        else:
#            RM_coeff = None;
#        b, _, _ =  bounds.NeronSeveriBound.verify_curve(g = g, factorsRR_geom = factorsRR_geom, RM_coeff = RM_coeff);
#        return b;

#  TODO: These methods need updating: use single Magma functions
    def verify_saturated(self):
        self._sat_test_, self._sat_cert_ =  magma.VerifySaturated(self._geo_rep_list_, self._P_, nvals = 2)
        return self._sat_test_

    def verify_representation(self):
        self._rep_test_ = True
        for gen in self._geo_rep_dict_:
            genTan = gen['tangent']
            test, corresp = self.correspondence(genTan)
            if not test:
                self._rep_test_ = False
            else:
                gen['corresp'] = corresp
        return self._rep_test_

    def verify(self):
        return (self.verify_algebra() and self.verify_saturated() and self.verify_representation())

class Lattice:
    def __init__(self, Endo):
        self.X = Endo.X
        self.g = Endo.g
        self.F = Endo.F
        self._geo_rep_list_ = Endo._geo_rep_list_
        self._list_, self._sthash_ = magma.EndomorphismLattice(self._geo_rep_list_, nvals = 2)
        self._desc_ = desc_lattice(self._list_)
        self._sthash_ = desc_sthash(self._sthash_)

    def __repr__(self):
        return repr_lattice(self)

    def _calculate_dictionary_(self):
        if not hasattr(self, "_dict_"):
            self._dict_ = dict_lattice(self._list_)

    def full(self):
        self._calculate_dictionary_()
        return self._dict_

    def pretty_print(self):
        return pretty_print_lattice_description(self._desc_, self.g)

#  TODO: These methods need updating: use single Magma functions
class Decomposition:
    def __init__(self, Endo):
        self.X = Endo.X
        self.g = Endo.g
        self.F = Endo.F
        self._P_ = Endo._P_
        idems, self.field = magma.IdempotentsFromLattice(Endo._lat_._list_, nvals = 2)
        self._facs_ = [ ]
        for idem in idems:
            fac = dict()
            fac['field'] = self.field
            fac['idem'] = dict_gen(idem)
            lat, proj = magma.ProjectionFromIdempotent(self._P_, idem, nvals = 2)
            fac['proj'] = dict_gen(proj);
            fac['factor'] = { 'analytic': lat }
            self._facs_.append(fac)

    def __repr__(self):
        return repr_decomposition(self)

    def full(self):
        return self._facs_

    def _calculate_factors_(self):
        for fac in self._facs_:
            if not 'algebraic' in fac['factor'].keys():
                fac['factor']['algebraic'] = magma.FactorReconstruct(self._P_, fac['factor']['analytic'], fac['proj']['approx'], fac['proj']['homology'], fac['field'])

    def factors(self):
        self._calculate_factors_()
        return [ fac['factor']['algebraic'] for fac in self._facs_ ]

    def _factors_desc_(self):
        self._calculate_factors_()
        return [ sagify_description(magma.FactorDescription(fac, self.F)) for fac in self.factors() ]
