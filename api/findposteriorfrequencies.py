import numpy as np
import pandas as pd

from api.calculate_frequencies import Frequencies

def findposteriorfrequencies(x: int, tempdata: np.ndarray, maxMOI: int, frequencies_RR, rand: np.random.RandomState):
    nalleles = frequencies_RR.lengths[x]

    # hard coded table() function from R
    data = tempdata[:, x*maxMOI: (x+1)*maxMOI].astype(int)
    data_1d_array = data.flatten()
    data_1d_array = data_1d_array[data_1d_array != 0]
    data_1d_array = data_1d_array[data_1d_array <= nalleles]

    data_unique, data_counts = np.unique(data_1d_array, return_counts=True)
    # Start as ones for the frequency prior
    counts_table = np.ones(nalleles)
    counts_table[data_unique - 1] += data_counts

    frequencies_RR.matrix[x, :nalleles] = rand.dirichlet(counts_table, 1)
