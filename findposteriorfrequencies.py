import numpy as np


def findposteriorfrequencies(x, tempdata, maxMOI, frequencies_RR):
    # data = tempdata[:, 1:maxMOI+(x)*maxMOI]
    data = tempdata[:,np.arange(1, maxMOI+1) + (x * maxMOI) - 1]
    nalleles = frequencies_RR[0][x]

    freq_prior_alpha = [1] * nalleles

    # hard coded table() function from R
    dictionary = [[i, 0] for i in range(1, nalleles + 1)]
    data_1d_array = data.flatten()
    data_1d_array = data_1d_array[data_1d_array != 0]

    for d in data_1d_array:
        if d <= nalleles:
            dictionary[d.astype(int) - 1][1] = dictionary[d.astype(int) - 1][1] + 1

    table = []
    for entry in dictionary:
        table.append(entry[1])

    freq_prior_alpha = np.asarray(freq_prior_alpha)
    table = np.asarray(table)

    freq_posterior_alpha = freq_prior_alpha + table
    freq_posterior_alpha = freq_posterior_alpha.tolist()
    frequencies_RR[1][x, 0:nalleles] = np.random.mtrand.dirichlet(freq_posterior_alpha, 1)
