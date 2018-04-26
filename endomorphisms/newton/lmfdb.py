#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

import pymongo

_C = None

def makeDBconnection():
    global _C

    _C = pymongo.MongoClient("m0.lmfdb.xyz");
    _C.admin.authenticate("lmfdb","lmfdb")


def getDBconnection():
    if _C is None:
        makeDBconnection()
    return _C

def get_curve(label):
    C = getDBconnection()
    return C.genus2_curves.curves.find({ u'label' : label}).next()

def get_endo(label):
    C = getDBconnection()
    return C.genus2_curves.endomorphisms.find({ u'label' : label}).next()
