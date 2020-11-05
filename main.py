import os

import numpy as np

from api.algorithm_instance import AlgorithmInstance, AlgorithmResults


# user inputs
inputfile = "Angola2017_example.xlsx"
locirepeats = np.array([2,2,3,3,3,3,3])
nruns = 10

# Import/set up data
test_run = AlgorithmInstance(inputfile, locirepeats)

# calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
record_interval = np.ceil(nruns / 1000)
burnin = np.ceil(nruns * 0.25)

results = test_run.run_algorithm(nruns, burnin, record_interval)

# Output summary data to .csv files
posterior_recrudescence_distribution_df, probability_of_recrudescence_df = results.get_summary_stats()
for site_name, posterior_df in results.run_posterior_dfs.items():
    posterior_df.to_csv(f'{site_name}_posterior.csv')
for site_name, summary_df in results.run_summary_stat_dfs.items():
    summary_df.to_csv(f'{site_name}_summary_statistics.csv')

posterior_recrudescence_distribution_df.to_csv('microsatellite_correction.csv', index=False)
probability_of_recrudescence_df.to_csv('probability_of_recrudescence.csv', index=False)
