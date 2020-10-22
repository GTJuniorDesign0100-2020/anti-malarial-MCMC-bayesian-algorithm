import os

import numpy as np

from algorithm_instance import AlgorithmInstance


# user inputs
inputfile = "Angola2017_example.xlsx"
locirepeats = np.array([2,2,3,3,3,3,3])
nruns = 1000

# Import/set up data
test_run = AlgorithmInstance(inputfile, locirepeats)

# calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
record_interval = np.ceil(nruns / 1000)
burnin = np.ceil(nruns * 0.25)

test_run.run_algorithm(nruns, burnin, record_interval)
