import numpy as np
import pandas as pd



# Todo:
#   Table layout can be handled by a configuration file/object.
#   Consider moving definitions related to the organization
#   of data frames into an object to decouple this code
#   from data organization.

def onLoad(genotypedata_latefailures,additional_genotypedata):
    site_names = siteNames()
    state_classification_all = np.array()    # dimensionality?
    state_parameters_all = np.array()        # dimensionality?
    ids_all = np.array()                    # 1-D vector

    for site in site_names:
        # jobname = site | unecessary declaration.

        # all rows from the table which have the current site.
        site_column = 1
        condition = (genotypedata_latefailures[:, site_column] == site)
        genoTypeData_RR = genotypedata_latefailures[condition]

        # new table, might have new row.
        site_column = 1
        condition = (additional_genotypedata[:, site_column] == site)
        additional_neutral = additional_genotypedata[condition]

        # Remember that R indexes from 1, not zero.
        if additional_neutral.shape[0] > 0:
            range = np.arange(1, additional_neutral.shape[0])
            str_range = np.array2string(range, separator='')
            str_range = str_range[1:len(str_range) - 1]
            str_append = 'Additional_' + str_range
            # This is where I am beginning to see an issue.
            # At this point in the R script, we set the Sample element's ID to 'str_append'.
            # That doesn't work here. Why not? Because we're working with a numpy array, not a DataFrame.
            # In short, numpy array wont HAVE elements, its an array, not an obj.

        # I THINK that mcmc.r is generating all this when it is run via source('mcmc.r') in the r script.
        (state_classification, state_parameters, ids) = runMCMC()

        # Not sure about concatenate vs append. Must be tested to confirm this concatenates rows of two similar tables.
        state_classification_all = np.append(state_classification, axis=0)
        state_parameters_all = np.append(state_parameters, axis=0)
        ids_all = np.append(ids)

    posterior_distribution_or_recrudescence = np.array()
    #(?) Set the name of the first (second?) column of posterior table to 'ID'

    # Takes the mean of each column for every row in a 2-dim array/matrix.
    # Behavior on empty tables needs to be evaluated.
    probability_of_recrudescence = np.mean(state_classification_all, axis=1)

    # Code says to generate a histogram at this point. Is that necessary? Doesn't seem like anything
    # Refrences that historgram or saves it for later viewing. If needed, use matplot package in the future.
    # Make histogram here.

    #Todo: Save to CSV files using Pandas.

    return

# place holder for a global variable. Use it to wrap the logic to access a vector of site names for now.
def siteNames():
    all_sites = np.array()
    return all_sites

# runs the MCMC script.
def runMCMC():
    state_classification = np.array()
    state_parameters = np.array()
    ids = np.array()
    return state_classification, state_parameters, ids

