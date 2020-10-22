import typing

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
        structures needed to run the algorithm
        '''
        input_file_parser.parseFile(input_file_path)
