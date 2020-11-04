from typing import IO, Union

class DataFileParser:
    '''
    Handles parsing the given file type and extracting data from it in a format
    usable by the MCMC algorithm
    '''

    @classmethod
    def parse_file(cls, input_file: Union[str, IO]):
        '''
        Attempts to parse the given file; returns a set of data structures
        describing the information in the file. If the file could not be parsed,
        throws an exception describing the problem.

        :param input_file: The string path OR file object for the data file to
        attempt parsing
        :return:
        '''
        # TODO: Implement this
