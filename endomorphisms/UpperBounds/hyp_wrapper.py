#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.

from utils import get_frob_list_HyperellipticCurve
from upper_bounds import endomorphisms_upper_bound, RR_upper_bound

def hyperelliptic_endomorphisms_upper_bound(f, h, bound = 100):
    frob_list = get_frob_list_HyperellipticCurve(f, h, bound);
    return endomorphisms_upper_bound(frob_list);

def hyperelliptic_RR_upper_bound(f, h, bound = 100):
    frob_list = get_frob_list_HyperellipticCurve(f, h, bound);
    return RR_upper_bound(frob_list);
