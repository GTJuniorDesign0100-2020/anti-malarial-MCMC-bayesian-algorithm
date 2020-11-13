import os
import sys

import numpy as np
import pandas as pd
import pickle
import pytest

'''
Saves the state of arguments from a switch_hidden call.
This allows a user to debug the final value between different
implementations.

Typically saves to the root level of the code base.
'''
def save_switch_state(x, nloci, maxMOI, alleles_definitions_RR, state, filename):
    packagedState = (x, nloci, maxMOI, alleles_definitions_RR, state);
    file = open(filename + str(".pickle"), 'wb')
    pickle.dump(packagedState, file)
    return

'''
Loads the state of arguments from a switch_hidden call.
This allows a user to debug the final value between different
implementations.

Loads from the local workspace, which at time of righting,
is the test folder. NOTE: That means that you'll need to copy the
file save at the root level into the test level to use it.
'''
def load_switch_state(filename):
    file = open(filename + str(".pickle"), 'rb')
    packagedState = pickle.load(file)
    return packagedState

def test_runTime(filename="debug_switch_state"):
    packaged_state = load_switch_state(filename);
    x = packaged_state[0]
    nloci = packaged_state[1]
    maxMOI = packaged_state[2]
    alleles_definitions_RR = packaged_state[3]
    state = packaged_state[4]



    return
