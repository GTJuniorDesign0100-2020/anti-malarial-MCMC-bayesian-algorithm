import os

import numpy as np

from api.algorithm_instance import AlgorithmInstance, AlgorithmResults

if __name__ == '__main__':
    # user inputs
    inputfile = "Angola2017_example.xlsx"
    locirepeats = np.array([2,2,3,3,3,3,3])
    nruns = 1000

    # Import/set up data
    test_run = AlgorithmInstance(inputfile, locirepeats)

    # calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
    record_interval = np.ceil(nruns / 1000)
    burnin = np.ceil(nruns * 0.25)

    results = test_run.run_algorithm(nruns, burnin, record_interval, is_verbose=True)

    # Output summary data to .csv files
    file_content = results.get_output_file_text()
    for filename, text in file_content.items():
        csv_file = open(filename, 'w')
        csv_file.write(text)
        csv_file.close()
