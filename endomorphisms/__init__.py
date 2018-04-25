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

from sage.all import magma

import os
import inspect
filename = inspect.getframeinfo(inspect.currentframe())[0];
__magmapath__ = os.path.dirname(filename) + "/magma/"
magma.AttachSpec(__magmapath__ + 'spec')

from Curves import mHyperellipticCurve, mPlaneCurve
from Wrapper import EndomorphismData, OverField, Lattice, Decomposition
from UpperBounds import endomorphisms_upper_bound, RR_upper_bound, hyperelliptic_endomorphisms_upper_bound, hyperelliptic_RR_upper_bound
