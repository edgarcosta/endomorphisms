"""
 *  Initialization
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

if not '__endodir__' in globals():
    raise ImportError("Please set a value for __endodir__")

import os
cur = os.getcwd()

os.chdir(__endodir__)
magma.chdir(__endodir__)

magma.load('~/.magmarc')
magma.AttachSpec('spec')

load('heuristic/Curves.sage')
load('heuristic/Relative.sage')
load('heuristic/Reprs.sage')
load('heuristic/Dicts.sage')
load('heuristic/PrettyPrint.sage')
load('Wrapper.sage')

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

os.chdir(cur)
magma.chdir(cur)
