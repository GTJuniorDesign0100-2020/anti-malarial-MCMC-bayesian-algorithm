import numpy as np
import pandas as pd

# Use this when encapsulating the code for this script.
# def onLoad(genotypedata_latefailures,additional_genotypedata):
#     site_names = siteNames()
#     state_classification_all = pd.DataFrame()
#     state_parameters_all = pd.DataFrame()
#     ids_all = np.array()  # 1-D vector
#
#     for site in site_names:
#         # jobname = site | unecessary declaration.
#
#         # all rows from the table which have the current site.
#         site_column = 1
#         condition = (genotypedata_latefailures[:, site_column] == site)
#         genoTypeData_RR = genotypedata_latefailures[condition]
#
#         # new table, might have new row.
#         site_column = 1
#         condition = (additional_genotypedata[:, site_column] == site)
#         additional_neutral = additional_genotypedata[condition]
#
#         # Remember that R indexes from 1, not zero.
#         if additional_neutral.shape[0] > 0:
#             range = np.arange(1, additional_neutral.shape[0])
#             str_range = np.array2string(range, separator='')
#             str_range = str_range[1:len(str_range) - 1]
#             str_append = 'Additional_' + str_range
#             # Big Problem:
#             # At this point in the R script, we set the Sample element's ID to 'str_append'.
#             # That doesn't work here. Why not? Because we're working with a numpy array, not a DataFrame.
#             # In short, numpy array wont HAVE elements, its an array, not an obj.
#
#         # I THINK that mcmc.r is generating all this when it is run via source('mcmc.r') in the r script.
#         (state_classification, state_parameters, ids) = runMCMC()
#
#         # Not sure about concatenate vs append. Must be tested to confirm this concatenates rows of two similar tables.
#         state_classification_all = state_classification_all.append(state_classification, ignore_index=True)
#         state_parameters_all = state_parameters_all.append(state_parameters, ignore_index=True)
#         ids_all = np.append(ids)
#
#     # Okay, we are now going to put together the PDOR table by row binding/column binding our previous tables.
#
#     # In the original script, this is a cbind2 then a rbind.
#     # I think this simpler translation accomplishes the same effect.
#     posterior_distribution_of_recrudescence = pd.DataFrame.concat([ids_all, state_classification_all])
#
#     # Set the name of the first column of PDOR table to 'ID'
#     # I don't think that's necessary, since everything gets an ID as is in Panda DataFrames,
#     # but we can add if necessary.
#
#     # Will this operation work? I read that numpy operations work over DataFrames,
#     # but I need to check this behavior to confirm.
#     probability_of_recrudescence = np.mean(state_classification_all, axis=1)
#
#     # Code says to generate a histogram at this point. Is that necessary? Doesn't seem like anything
#     # Refrences that historgram or saves it for later viewing. If needed, use matplot package in the future.
#     # Make histogram here.
#
#     # Save csv's
#     posterior_distribution_of_recrudescence.write_csv('microsatellite_correction.csv', index=False)
#     probability_of_recrudescence.write_csv( "probability_of_recrudescence.csv", index=False)
#     return

# place holder for a global variable. Use it to wrap the logic to access a vector of site names for now.
def siteNames():
    all_sites = np.array()
    return all_sites

# runs the MCMC script.
# are these numpy arrays, or data frames?
def runMCMC():
    state_classification = pd.DataFrame()
    state_parameters = pd.DataFrame()
    ids = np.array()
    return state_classification, state_parameters, ids

def runMCMCscript():
    import mcmc.py
    return

####################################
#   Soup Style Start.
#   Don't do this long term lol.
####################################
site_names = siteNames()
state_classification_all = pd.DataFrame()
state_parameters_all = pd.DataFrame()
ids_all = np.array()  # 1-D vector

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
        # Big Problem:
        # At this point in the R script, we set the Sample element's ID to 'str_append'.
        # That doesn't work here. Why not? Because we're working with a numpy array, not a DataFrame.
        # In short, numpy array wont HAVE elements, its an array, not an obj.

    # I THINK that mcmc.r is generating all this when it is run via source('mcmc.r') in the r script.
    runMCMCscript()

    # Not sure about concatenate vs append. Must be tested to confirm this concatenates rows of two similar tables.
    state_classification_all = state_classification_all.append(state_classification, ignore_index=True)
    state_parameters_all = state_parameters_all.append(state_parameters, ignore_index=True)
    ids_all = np.append(ids)

# Okay, we are now going to put together the PDOR table by row binding/column binding our previous tables.

# In the original script, this is a cbind2 then a rbind.
# I think this simpler translation accomplishes the same effect.
posterior_distribution_of_recrudescence = pd.DataFrame.concat([ids_all, state_classification_all])

# Set the name of the first column of PDOR table to 'ID'
# I don't think that's necessary, since everything gets an ID as is in Panda DataFrames,
# but we can add if necessary.

# Will this operation work? I read that numpy operations work over DataFrames,
# but I need to check this behavior to confirm.
probability_of_recrudescence = np.mean(state_classification_all, axis=1)

# Code says to generate a histogram at this point. Is that necessary? Doesn't seem like anything
# Refrences that historgram or saves it for later viewing. If needed, use matplot package in the future.
# Make histogram here.

# Save csv's
posterior_distribution_of_recrudescence.write_csv('microsatellite_correction.csv', index=False)
probability_of_recrudescence.write_csv( "probability_of_recrudescence.csv", index=False)

####################################
#   Soup Style Script End.
#   Don't do this long term lol.
####################################