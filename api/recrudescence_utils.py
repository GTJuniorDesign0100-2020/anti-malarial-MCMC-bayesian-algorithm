'''
Contains general utility functions used for working with recrudscence data in
the malaria recrudescence algorithm

NOTE: When possible, many of these functions should be refactored into a more
specific module to avoid throwing everything in here
'''

import numpy as np
import pandas as pd


def get_sample_ids(genotypedata: pd.DataFrame, search_text: str) -> np.ndarray:
    '''
    Returns a numpy array of all the unique "Sample ID"s in the dataframe
    whose name contains the matching text

    NOTE: Currently assumes space in ID name
    TODO: Could potentially be refactored to a more general function?

    :param genotypedata: The data to search
    :param search_text: The substring to look for in the
    '''
    matching_sample_names = genotypedata[
            genotypedata['Sample ID'].str.contains(search_text)
        ]['Sample ID']
    # Remove day from sample
    sample_id_names = matching_sample_names.str.split(' ', n=1).str.get(0)
    return pd.unique(sample_id_names)


def nan_array(shape: tuple) -> np.ndarray:
    '''
    Returns a numpy array of the given size filled entirely with NaN values
    :param size: A tuple specifying the desired shape for the numpy array
    :return: A numpy array with the given shape filled with np.nan
    '''
    return np.full_like(np.empty(shape=shape), np.nan)
