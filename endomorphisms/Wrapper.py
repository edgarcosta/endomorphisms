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
            self._P_ = magma.PeriodMatrix(self.X)

    def __repr__(self):
        return repr_endomorphism_data(self)

    def lattice(self):
        if not hasattr(self, "_lat_list_"):
            self._lat_list_ = magma.HeuristicEndomorphismLattice(self.X)
            self._lat_desc_ = desc_lattice(self._lat_list_)
            self._lat_dict_ = dict_lattice(self._lat_list_)
        return pretty_print_lattice_description(self._lat_desc_, self.g)

    def decomposition(self):
        if not hasattr(self, "_dec_"):
            self._dec_ = magma.HeuristicJacobianFactors(self.X)
        return self._dec_

    def verify_upper_bound(self):
        if not hasattr(self, "_test_upper_"):
            self._test_upper_ = True
        return self._test_upper_

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

    def verify_lower_bound(self):
        if not hasattr(self, "_test_lower_"):
            self._test_lower_, self._test_lower_cert_ = magma.CertifiedEndomorphismAlgebra(self.X, Geometric = True, nvals = 2)
        return self._test_lower_

    def verify(self):
        test_upper = self.verify_upper_bound()
        test_lower = self.verify_lower_bound()
        return test_upper and test_lower
