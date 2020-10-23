'''
Contains general utility functions used for working with recrudscence data in
the malaria recrudescence algorithm

NOTE: When possible, many of these functions should be refactored into a more
specific module to avoid throwing everything in here
'''

import pandas as pd


def get_sample_ids(genotypedata: pd.DataFrame, search_text: str):
    '''
    Returns a numpy array of all the unique "Sample ID"s in the dataframe
    whose name contains the matching text

    TODO: Currently assumes underscore remains in ID name
    TODO: Could potentially be refactored to a more general function?

    :param genotypedata: The data to search
    :param search_text: The substring to look for in the
    '''
    matching_sample_names = genotypedata[
            genotypedata['Sample ID'].str.contains(search_text)
        ]['Sample ID']
    # Remove day from sample
    sample_id_names = matching_sample_names.str.rsplit('_', n=1).str.get(0)
    return pd.unique(sample_id_names)
