import numpy as np
import pandas as pd

from api.data_file_parser import DataFileParser
import api.recrudescence_utils as recrudescence_utils

class RecrudescenceFileParser(DataFileParser):
    '''
    Parses the given malaria test data file and returns the data. The data
    is returned in 2 parts:
    -   A set of genotypedata containing a set of initial samples and a day of
        failure when those samples returned positive for malaria again
    -   A set of background "additional" geotypedata composed of only single
        samples, and used for determining allele frequency, etc., sametimes
        known as the "background day 0 samples"
    '''

    @classmethod
    def parse_file(cls, input_file_path: str):
        '''
        Attempts to parse the given malaria test data file and returns data
        structures usable by the MCMC algorithm; if this fails, throws an
        exception describing what went wrong (and where, if the file itself
        was formatted incorrectly)

        :param input_file_path: The string path to the Excel data file to
        attempt parsing
        :return: Tuple of (genotypedata_latefailures, additional_genotypedata),
        with the former being the data for day 0/day of failure samples and the
        latter being background samples
        # TODO: Add proper exceptions/error handling
        '''
        raw_file_info = pd.ExcelFile(input_file_path)
        genotypedata_latefailures = cls._get_genotype_late_failures(raw_file_info)
        additional_genotypedata = cls._get_additional_genotype_data(raw_file_info)

        return genotypedata_latefailures, additional_genotypedata

    @classmethod
    def _get_genotype_late_failures(cls, raw_file_info: pd.ExcelFile):
        '''
        Returns a dataframe with all the sample data from the "Late Treatment
        Failures" tab of the Excel data spreadsheet

        :param raw_file_info: The input data file, opened as an Excel sheet in
        pandas
        :return: TODO: Describe structure
        '''
        genotypedata_latefailures = raw_file_info.parse(
            'Late Treatment Failures', skiprows=3)
        cls._replace_missing_data(genotypedata_latefailures)

        # recode sample names so that each pair has a " Day 0" and a " Day Failure"
        genotypedata_latefailures['Sample ID'] = genotypedata_latefailures['Sample ID'].str.replace('_D0$', ' Day 0')
        genotypedata_latefailures['Sample ID'] = genotypedata_latefailures['Sample ID'].str.replace('_D[0-9]+$', ' Day Failure')

        # verify each sample has a Day 0 and a Day of Failure
        day_0_ids = recrudescence_utils.get_sample_ids(genotypedata_latefailures, 'Day 0')
        day_fail_ids = recrudescence_utils.get_sample_ids(genotypedata_latefailures, 'Day Failure')

        if day_0_ids.size != day_fail_ids.size:
            # TODO: Use more specific exception
            raise Exception('Error - each sample must have day 0 and day of failure data')

        return genotypedata_latefailures

    @classmethod
    def _get_additional_genotype_data(cls, raw_file_info: pd.ExcelFile):
        '''
        Returns a dataframe with all the sample data from the "Additional" tab
        of the Excel data spreadsheet

        :param raw_file_info: The input data file, opened as an Excel sheet in
        pandas
        :return: TODO: Describe structure
        '''
        additional_genotypedata = raw_file_info.parse('Additional', skiprows=3)
        cls._replace_missing_data(additional_genotypedata)

        # recode sample names to " Day 0" and " Day Failure"
        additional_genotypedata['Sample ID'] = additional_genotypedata['Sample ID'].str.replace('_D0$', ' Day 0')
        additional_genotypedata['Sample ID'] = additional_genotypedata['Sample ID'].str.replace('_D[0-9]+$', ' Day Failure')

        # TODO: Implement this!
        return additional_genotypedata

    @classmethod
    def _replace_missing_data(cls, df: pd.DataFrame):
        '''
        Replaces any missing data in the dataframe with NaN values (i.e.
        numpy.nan)

        :param df: The pandas dataframe to replace the values in
        :return: The dataframe with the missing values replaced (although it
        should update the original dataframe anyway)
        '''
        df.replace(0, np.nan, inplace=True)
        df.replace('0', np.nan, inplace=True)
        df.replace('N/A', np.nan, inplace=True)
        df.replace('NA', np.nan, inplace=True)
        df.replace('-', np.nan, inplace=True)
        return df
