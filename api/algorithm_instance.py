from collections import OrderedDict
from typing import IO, List, Union

import numpy as np
import pandas as pd

from api.algorithm_site_instance import AlgorithmSiteInstance, SavedState
from api.data_file_parser import DataFileParser
from api.recrudescence_file_parser import RecrudescenceFileParser


class AlgorithmResults:
    '''
    Holds the final results of running the malaria recrudescence algorithm
    '''

    def __init__(self, site_names):
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
        self.saved_classification_all = {}
        self.saved_parameters_all = {}
        self.sample_ids = np.array([])
        self.site_sample_ids = {}
        self.run_posterior_dfs = {}
        self.run_summary_stat_dfs = {}
        self.site_names = site_names

    def update_results(self, site_results: SavedState, site_name: str):
        '''
        Add the site's results to the overall algorithm results
        '''
        # TODO: Determine output array is correct shape
        self.saved_classification_all[site_name] = pd.DataFrame(site_results.classification)
        self.saved_parameters_all[site_name] = pd.DataFrame(site_results.parameters)
        self.sample_ids = np.append(self.sample_ids, site_results.ids)

        self.site_sample_ids[site_name] =  site_results.ids
        self.run_posterior_dfs[site_name] = pd.DataFrame(site_results.posterior_df)
        self.run_summary_stat_dfs[site_name] = pd.DataFrame(site_results.summary_stats_df)

    def get_summary_stats(self):
        '''
        Calculate the final recrudescence probabilities and distributions for
        each sample in the dataset, across all sites
        :return: A tuple of 2 pandas dataframes,
        (posterior_recrudescence_distribution, probability_of_recrudescence)
        '''
        posterior_recrudescence_distribution = pd.DataFrame()
        ids = pd.DataFrame()
        for site in self.site_names:
            posterior_recrudescence_distribution = posterior_recrudescence_distribution.append(self.saved_classification_all[site])
            ids = ids.append(pd.DataFrame(self.site_sample_ids[site]))

        data = posterior_recrudescence_distribution.copy()
        new_col = []
        for col in posterior_recrudescence_distribution.columns:
            if type(col) is int:
                new_col.append("V" + str(col))
            else:
                new_col.append(col)
        posterior_recrudescence_distribution.columns = new_col
        posterior_recrudescence_distribution.insert(0, "Sample ID", ids, True)
        # Get the mean of each row
        probability_of_recrudescence = pd.DataFrame(np.mean(data, axis=1), columns = ["Recrudescence Probability"])
        probability_of_recrudescence.insert(0, "Sample ID", ids, True)

        return posterior_recrudescence_distribution, probability_of_recrudescence

    def get_output_file_text(self) -> dict:
        '''
        Gets the text of each output .csv file and stores it in a dictionary,
        with the filename as the key
        :return: A dictionary with the filename of each .csv file as the key
        and the content of the file as its value
        '''
        output_files = {}
        extended_ids = {key: [] for key in self.site_sample_ids.keys()}
        for site in self.site_sample_ids:
            for id in self.site_sample_ids[site]:
                extended_ids[site].append(str(id) + " Day 0")
                extended_ids[site].append(str(id) + " Day Failure")

        for site_name, posterior_df in self.run_posterior_dfs.items():
            filename = f'{site_name}_posterior.csv'
            posterior_df.insert(0, "Sample ID", extended_ids[site_name], True)
            output_files[filename] = posterior_df.to_csv(index=False, line_terminator='\n')

        desc = []
        loci = []
        for site_name, summary_df in self.run_summary_stat_dfs.items():
            if not desc:
                loci = summary_df.index.values[2:-1]
                num_loci = int(len(loci)/2)
                desc.append("Probability of Missing Allele")
                desc.append("Error Rate")
                for i in range(num_loci):
                    desc.append("Frequency of most Common Allele")
                for i in range(num_loci):
                    desc.append("Simpson's Diversity Index")
                desc.append("Average Simpson's Index")
            filename = f'{site_name}_summary_statistics.csv'
            summary_df.insert(0, "Description", desc, True)
            locus_names = list(summary_df.index.values[0:-1])
            locus_names.append(" ")
            summary_df.insert(1, "Locus Name", locus_names, True)
            summary_df.columns = ["Description", "Locus Name", "Mean with Interquartile range (25%, 75%)"]
            output_files[filename] = summary_df.to_csv(index=False, line_terminator='\n')
        for site_name, param_df in self.saved_parameters_all.items():
            filename = f'{site_name}_state_parameters.csv'
            if (len(desc) % 2 != 0):
                desc = desc[0:-1]
            param_df.insert(0, "Description", desc, True)
            loci1 = ['q', 'd']
            res = [y for x in [loci1, loci] for y in x]
            param_df.insert(1, "Locus Name", res, True)
            new_col = []
            for col in param_df.columns:
                if type(col) is int:
                    new_col.append("V" + str(col))
                else:
                    new_col.append(col)
            param_df.columns = new_col
            output_files[filename] = param_df.to_csv(index=False, line_terminator='\n')

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df = self.get_summary_stats()
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
        self.site_names = pd.unique(genotypedata['Site'])
        self.algorithm_instances = []
        for site_name in self.site_names:
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
        overall_results = AlgorithmResults(self.site_names)

        for site_name, algo_instance in self.algorithm_instances:
            site_result = algo_instance.run_algorithm(
                site_name,
                nruns,
                burnin,
                record_interval,
                seed)
            # save site results
            overall_results.update_results(site_result, site_name)

        return overall_results
