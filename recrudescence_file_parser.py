import numpy as np
import pandas as pd

from data_file_parser import DataFileParser


class RecrudescenceFileParser(DataFileParser):
    '''
    Parses the given malaria test data file and returns
    '''

    @classmethod
    def parseFile(cls, input_file_path: str):
        '''
        Attempts to parse the given malaria test data file and returns data
        structures usable by the MCMC algorithm; if this fails, throws an
        exception describing what went wrong (and where, if the file itself
        was formatted incorrectly)

        :param input_file_path: The string path to the Excel data file to
        attempt parsing
        :return: Tuple of (genotypedata_latefailures, additional_genotypedata)
        # TODO: Figure out what these actually mean?
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
        # TODO: Implement this!
        return None

    @classmethod
    def _get_additional_genotype_data(cls, raw_file_info: pd.ExcelFile):
        '''
        Returns a dataframe with all the sample data from the "Late Treatment
        Failures" tab of the Excel data spreadsheet

        :param raw_file_info: The input data file, opened as an Excel sheet in
        pandas
        :return: TODO: Describe structure
        '''
        # TODO: Implement this!
        return None

    @classmethod
    def _replace_missing_data(cls, dataframe: pd.DataFrame):
        '''
        Replaces any missing data in the dataframe with NaN values (i.e.
        numpy.nan)

        :param dataframe: The pandas dataframe to replace the values in
        :return: The dataframe with the missing values replaced (although it
        should update the original dataframe anyway)
        '''
        dataframe.replace(0, np.nan)
        dataframe.replace('0', np.nan)
        dataframe.replace('N/A', np.nan)
        dataframe.replace('NA', np.nan)
        dataframe.replace('-', np.nan)
        return dataframe
