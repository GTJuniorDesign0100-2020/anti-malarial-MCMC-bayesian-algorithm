import typing

import numpy as np
import pandas as pd

from data_file_parser import DataFileParser
from recrudescence_file_parser import RecrudescenceFileParser


class AlgorithmInstance:
    '''
    Handles setting up and running an instance of the MCMC recrudescence
    algorithm on the given malaria test data file
    '''

    def __init__(
        self,
        input_file_path: str,
        locirepeats: typing.List[int],
        input_file_parser: DataFileParser=RecrudescenceFileParser):
        '''
        Parses the given malaria test data file and sets up the initial data
        structures needed to run the algorithm for each "arm"/site location in
        the file
        :param input_file_path: The string path to the data file
        :param locirepeats: TODO:
        :param input_file_parser: (optional) The parser specifying how to read/
        interpret the data
        '''
        genotypedata_latefailures, additional_genotypedata = input_file_parser.parse_file(input_file_path)

        site_names = pd.unique(genotypedata_latefailures['Site'])
        self.algorithm_instances = []
        for site_name in site_names:
            site_genotypedata_RR = genotypedata_latefailures[
                    genotypedata_latefailures['Site'] == site_name
                ].drop(columns='Site')
            site_additional_neutral = additional_genotypedata[
                    additional_genotypedata['Site'] == site_name
                ].drop(columns='Site')

            self._replace_sample_names(site_additional_neutral, 'Additional_')

            site_instance = AlgorithmSiteInstance(
                site_genotypedata_RR,
                site_additional_neutral,
                locirepeats)
            self.algorithm_instances.append((site_name, site_instance))

    @classmethod
    def _replace_sample_names(cls, samples_df: pd.DataFrame, text: str='Additional_'):
        '''
        Replaces all of the "Sample ID"s in the given dataframe with their
        index, preceded by a piece of text (e.g. the 1st sample with the
        default text would be named "Additional_0")

        :param samples_df: The dataframe to replace the sample names in
        :param text: (optional) The text to prepend before the sample's index
        :return: The dataframe with the replaced sample ID names (although it
        should modify the original samples_df)
        '''
        new_sample_names = np.array(
            [f'{text}{i}' for i in range(samples_df.shape[0])])
        samples_df['Sample ID'] = new_sample_names
        return samples_df

    def run_algorithm(self, nruns: int=1000, burnin: int=100, record_interval: int=10, seed=None):
        '''
        Runs the actual MCMC algorithm for each site, and returns the combined
        data for all of the sites
        '''
        # Set up data structures to hold results
        state_classification_all = pd.DataFrame()
        state_parameters_all = pd.DataFrame()
        ids_all = np.array([])
        # TODO: Finish implementing this!


class AlgorithmSiteInstance:
    '''
    Handles running the MCMC malaria recrudescence algorithm for a single
    "arm"/site's data
    '''

    def __init__(
        self,
        genotypedata_RR: pd.DataFrame,
        additional_neutral: pd.DataFrame,
        locirepeats: typing.List[int]):
        '''
        Sets up the initial data structures needed before running the algorithm
        '''
        # TODO: Implement this!
