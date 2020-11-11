from collections import OrderedDict, namedtuple
import multiprocessing
from typing import IO, List, Tuple, Union

import numpy as np
import pandas as pd

from api.algorithm_site_instance import AlgorithmSiteInstance, SavedState
from api.data_file_parser import DataFileParser
from api.recrudescence_file_parser import RecrudescenceFileParser


RunInputArgs = namedtuple('RunInputArgs', 'nruns burnin record_interval seed')


class AlgorithmResults:
    '''
    Holds the final results of running the malaria recrudescence algorithm
    '''

    def __init__(self):
        '''
        Set up data structures the hold the results, as follows:
        saved_classification_all - Holds the recrudescence/reinfection
        classifications given to each sample
        saved_parameters_all - Holds the hidden parameter values of each site as
        the algorithm ran
        ids_all - Lists the IDs of each sample in the entire dataset
        site_sample_ids - Lists the IDs of each sample in a given site
        run_posterior_dfs - For each site, holds the computed posterior
        distributions (held in pandas dataframes)
        run_summary_stats - For each site, holds the computed summary statistics
        (in a pandas dataframe)
        '''
        self.saved_classification_all = pd.DataFrame()
        self.saved_parameters_all = pd.DataFrame()
        self.sample_ids = np.array([])
        self.site_sample_ids = OrderedDict()
        self.run_posterior_dfs = {}
        self.run_summary_stat_dfs = {}

    def update_results(self, site_results: SavedState, site_name: str):
        '''
        Add the site's results to the overall algorithm results
        '''
        # TODO: Determine output array is correct shape
        self.saved_classification_all = self.saved_classification_all.append(
            pd.DataFrame(site_results.classification), ignore_index=True)
        self.saved_parameters_all = self.saved_parameters_all.append(
            pd.DataFrame(site_results.parameters), ignore_index=True)
        self.sample_ids = np.append(self.sample_ids, site_results.ids)

        self.site_sample_ids[site_name] = site_results.ids
        self.run_posterior_dfs[site_name] = site_results.posterior_df
        self.run_summary_stat_dfs[site_name] = site_results.summary_stats_df

    def get_summary_stats(self):
        '''
        Calculate the final recrudescence probabilities and distributions for
        each sample in the dataset, across all sites
        :return: A tuple of 2 pandas dataframes,
        (posterior_recrudescence_distribution, probability_of_recrudescence)
        '''
        posterior_recrudescence_distribution = pd.concat([
            pd.DataFrame(self.sample_ids), self.saved_classification_all])
        posterior_recrudescence_distribution.rename(columns={
            posterior_recrudescence_distribution.columns[0]: 'ID'
        }, inplace=True)
        # Get the mean of each row
        probability_of_recrudescence = np.mean(self.saved_classification_all, axis=1)

        return posterior_recrudescence_distribution, probability_of_recrudescence

    def get_output_file_text(self) -> dict:
        '''
        Gets the text of each output .csv file and stores it in a dictionary,
        with the filename as the key
        :return: A dictionary with the filename of each .csv file as the key
        and the content of the file as its value
        '''
        output_files = {}

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df = self.get_summary_stats()
        for site_name, posterior_df in self.run_posterior_dfs.items():
            filename = f'{site_name}_posterior.csv'
            output_files[filename] = posterior_df.to_csv(line_terminator='\n')
        for site_name, summary_df in self.run_summary_stat_dfs.items():
            filename = f'{site_name}_summary_statistics.csv'
            output_files[filename] = summary_df.to_csv(line_terminator='\n')

        output_files['microsatellite_correction.csv'] = posterior_recrudescence_distribution_df.to_csv(index=False, line_terminator='\n')
        output_files['probability_of_recrudescence.csv'] = probability_of_recrudescence_df.to_csv(index=False, line_terminator='\n')

        return output_files


class AlgorithmInstance:
    '''
    Handles setting up and running an instance of the MCMC recrudescence
    algorithm on the given malaria test data file (for all "arms"/sites in the
    data)
    '''

    def __init__(
        self,
        input_file: Union[str, IO],
        locirepeats: List[int],
        input_file_parser: DataFileParser=RecrudescenceFileParser):
        '''
        Parses the given malaria test data file and sets up the initial data
        structures needed to run the algorithm for each "arm"/site location in
        the file

        :param input_file: The string path OR file object for the data file to
        process
        :param locirepeats: TODO:
        :param input_file_parser: (optional) The parser specifying how to read/
        interpret the data
        '''
        genotypedata_latefailures, additional_genotypedata = input_file_parser.parse_file(input_file)

        self._initialize_site_instances(genotypedata_latefailures, additional_genotypedata, locirepeats)

    def _initialize_site_instances(self, genotypedata: pd.DataFrame, additional: pd.DataFrame, locirepeats: List[int]):
        '''
        Initializes the algorithm instances for each site in the processed
        data file
        :param locirepeats: TODO:
        '''
        site_names = pd.unique(genotypedata['Site'])
        self.algorithm_instances = []
        for site_name in site_names:
            # NOTE: "RR" stands for "recrudescence and/or reinfection"; it marks
            # datasets that deals specifically with day 0/day of failure info,
            # as opposed to background data
            site_genotypedata_RR = self._get_samples_from_site(
                genotypedata, site_name)
            site_additional_neutral = self._get_samples_from_site(
                additional, site_name)

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

    def run_algorithm(self, nruns: int=1000, burnin: int=100, record_interval: int=10, seed=None) -> AlgorithmResults:
        '''
        Runs the actual MCMC algorithm for each site, and returns the combined
        data for all of the sites
        TODO: Elaborate
        '''
        overall_results = AlgorithmResults()

        num_sites = len(self.algorithm_instances)
        with multiprocessing.Pool(processes=num_sites) as pool:
            run_inputs = [RunInputArgs(nruns, burnin, record_interval, seed)] * num_sites
            all_site_arguments = zip(self.algorithm_instances, run_inputs)
            site_results = pool.starmap(
                self._run_site_algorithm,
                all_site_arguments)

            for site_name, site_result in site_results:
                # save site results
                overall_results.update_results(site_result, site_name)

        return overall_results

    @classmethod
    def _run_site_algorithm(cls,
        site_instance: Tuple[str, AlgorithmSiteInstance],
        run_inputs: RunInputArgs) -> Tuple[str, SavedState]:
        '''
        Runs the algorithm for a single site instance, and returns the
        algorithm results along with the site name
        :param site_instance: The site to start running the algorithm for
        :run_inputs: The input arguments to the algorithm run
        :return: A tuple containing the site name and site algorithm results
        '''
        site_name, algo_instance = site_instance
        site_result = algo_instance.run_algorithm(
            site_name,
            run_inputs.nruns,
            run_inputs.burnin,
            run_inputs.record_interval,
            run_inputs.seed)
        return site_name, site_result
