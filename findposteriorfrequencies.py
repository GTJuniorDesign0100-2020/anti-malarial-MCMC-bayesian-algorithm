import numpy as np
import run_all_arms;
import mcmc;


def findposteriorfrequencies(x, tempdata, maxMOI, frequencies_RR):
    data = tempdata[:,1:maxMOI+(x-1)*maxMOI];
    nalleles = frequencies_RR[[1]][x];
    freq_prior_alpha = [1]*nalleles;

    # hard coded table() function from R 
    table = [[i, 0] for i in range(1,nalleles)];
    for d in data:
    	table[d-1][1] = table[d-1][1]+1;
    freq_posterior_alpha = freq_prior_alpha + table;

    frequencies_RR[[2]][x,1:nalleles] = numpy.random.mtrand.dirichlet(freq_posterior_alpha, 1);
	
