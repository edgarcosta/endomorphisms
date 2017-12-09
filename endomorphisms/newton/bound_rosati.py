#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.
from sage.all import Matrix

def bound_rosati(alpha_geo):
    H = Matrix(4,4)
    H[0, 2] = H[1, 3] = -1
    H[2, 0] = H[3, 1] = 1
    return (alpha_geo * H * alpha_geo.transpose() * H.inverse()).trace()/2
