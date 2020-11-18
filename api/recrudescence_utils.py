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

def get_locinames(genotypedata: pd.DataFrame):
    '''
    Returns a dictionary that has the names of locus without any duplicates
    the list will contain the names with the order that was in the genotypedata (left to right)

    :param genotypedata: The dataframe that has all names of locus with duplicates
    '''
    col_names = list(genotypedata.columns.values)[1:]
    prev_pos = None
    locinames = {}
    lociname_index = 0
    lociname_end_index = 0
    for name in col_names:
        res = name.split("_")
        if prev_pos == None:
            prev_pos = res[0]
        elif prev_pos != res[0]:
            locinames[lociname_index] = (prev_pos, lociname_end_index)
            prev_pos = res[0]
            lociname_index += 1
            lociname_end_index += 1
        else:
            lociname_end_index += 1
            if (lociname_end_index == len(col_names)-1):
                locinames[lociname_index] = (prev_pos, lociname_end_index)
    return locinames

def get_RawAlleles(genotypedata: pd.DataFrame, current_column: int, last_column: int):
    '''
    Given the dataframe, returns a list of all allele values 
    of the Last Treatment Failures and Additional tables for selected loci name
    (only from the columns between the initial value of current_column and last_column)


    :param genotypedata: The data table to retrieve allele values
    :param current_column: The current column index of the genotypedata
    :param last_column: The last column index that this method should retrieve allele values
    '''
    raw_alleles = []
    while (current_column <= last_column):
        raw_alleles += genotypedata.iloc[:, current_column+1].tolist()
        current_column += 1
    raw_alleles = [loci for loci in raw_alleles if str(loci) != 'nan']
    return raw_alleles, current_column