#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.

# exposes some of the functionality mentioned in Section 7.3 and Section 7.4

from .upper_bounds import endomorphisms_upper_bound, RR_upper_bound
assert endomorphisms_upper_bound
assert RR_upper_bound
from .hyp_wrapper  import hyperelliptic_endomorphisms_upper_bound, hyperelliptic_RR_upper_bound
assert hyperelliptic_endomorphisms_upper_bound
assert hyperelliptic_RR_upper_bound


