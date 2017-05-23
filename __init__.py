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

from sage.all import load

# by using inspect, we can set __endodir__ to the right path, and use heuristic_endomorphisms as a python module
if not '__endodir__' in globals():
    import os
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe())[0];
    __endodir__ = os.path.dirname(filename) + "/"

from sage.all import *
load(__endodir__ + "Initialize.sage");
