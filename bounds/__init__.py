"""
 *  Initialization of bound functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

from sage.all import load

import os
__boundsdir__ = os.getcwd() + '/bounds/'

from sage.all import *
magma.AttachSpec(__boundsdir__ + "spec")
load(__boundsdir__ + "constants.sage");
load(__boundsdir__ + "DiscriminantBound.sage")
load(__boundsdir__ + "TwistPolynomials.sage")
load(__boundsdir__ + "NonQM.sage");
load(__boundsdir__ + "GeometricallyIrreducible.sage");
load(__boundsdir__ + "EndomorphismRankBound.sage")
load(__boundsdir__ + "NonIsogenous.sage")
load(__boundsdir__ + "Genus2Factors.sage")
load(__boundsdir__ + "PointCounting.sage")
load(__boundsdir__ + "NonSquareCM.sage")
load(__boundsdir__ + "ProductsEC.sage")
