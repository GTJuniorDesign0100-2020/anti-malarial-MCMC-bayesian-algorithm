import typing

import numpy as np
import pandas as pd

from algorithm_site_instance import AlgorithmSiteInstance
from data_file_parser import DataFileParser
from recrudescence_file_parser import RecrudescenceFileParser


class AlgorithmInstance:
    '''
    Handles setting up and running an instance of the MCMC recrudescence
    algorithm on the given malaria test data file (for all "arms"/sites in the
    data)
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
            site_genotypedata_RR = self._get_samples_from_site(
                genotypedata_latefailures, site_name)
            site_additional_neutral = self._get_samples_from_site(
                additional_genotypedata, site_name)

            self._replace_sample_names(site_additional_neutral, 'Additional_')

            site_instance = AlgorithmSiteInstance(
                site_genotypedata_RR,
                site_additional_neutral,
                locirepeats)
            self.algorithm_instances.append((site_name, site_instance))

    @classmethod
    def _get_samples_from_site(cls, samples_df: pd.DataFrame, site_name: str):
        '''
        Returns a dataframe that only contains samples from the given site, with
        the site information removed from the dataframe

        :param sample_df: The dataframe containing the sample data
        :param site_name: The name of the site we want samples for
        :return: A new dataframe containing only samples from "site_name" with
        the "Site" column dropped
        '''
        return samples_df[samples_df['Site'] == site_name].drop(columns='Site')

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
        TODO: Elaborate
        '''
        # Set up data structures to hold results
        saved_classification_all = pd.DataFrame()
        saved_parameters_all = pd.DataFrame()
        ids_all = np.array([])

        for site_name, algo_instance in self.algorithm_instances:
            saved_classification, saved_params, ids = algo_instance.run_algorithm(
                site_name,
                nruns,
                burnin,
                record_interval,
                seed)

            # save site results
            saved_classification_all = pd.concat(saved_classification_all, saved_classification)
            saved_parameters_all = pd.concat(saved_parameters_all, saved_params)
            ids_all = np.append(ids_all, ids)

        return self._get_summary_stats(saved_classification_all, saved_parameters_all, ids_all)

    def _get_summary_stats(self, saved_classification, saved_parameters, ids):
        '''
        Calculate the final recrudescence probabilities and distributions for
        each sample
        TODO: Elaborate
        '''
        posterior_recrudescence_distribution = pd.concat(ids, saved_classification)
        posterior_recrudescence_distribution.rename(columns={
            posterior_recrudescence_distribution.columns[0]: 'ID'
        }, inplace=True)
        # Get the mean of each row
        probability_of_recrudescence = np.mean(saved_classification, axis=1)

        # NOTE: Not saving to .csv yet as in original code (make that a separate module)
        return posterior_recrudescence_distribution, probability_of_recrudescence
