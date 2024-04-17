import feelpp
from ._wave import *
import json

_waves = {
    'wave(2,1)': Wave2DP1,
    'wave(2,2)': Wave2DP2,
}


def get(dim=2, order=1, worldComm=None):
    """create a Wave operator

    """
    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()
    key = 'wave('+str(dim)+','+str(order)+')'
    if key not in _waves:
        raise RuntimeError('Wave'+key+' is not available')
    return _waves[key]()

def loadSpecs(jsonfile):
    # Reading the JSON file
    with open(jsonfile, 'r') as file:
        data = json.load(file)
    return data

