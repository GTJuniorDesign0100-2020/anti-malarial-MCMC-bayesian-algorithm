import os
import numpy as np
import Import_Microsatellite_Data
import run_all_arms

# set working directory
# os.chdir("\\\\cdc.gov\\private\\M326\\wif7\\Microsatellite\\GATech\\MRE");

#import gtools;


# user inputs
inputfile = "Angola2017_example.xlsx";
locirepeats = np.array([[2,2,3,3,3,3,3]]);
nruns = 1000;

# call script to import data
genotypedata_latefailures, additional_genotypedata = Import_Microsatellite_Data.onload()

# calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
record_interval = np.ceil(nruns / 1000);
burnin = np.ceil(nruns * 0.25);
run_all_arms.onload(genotypedata_latefailures, additional_genotypedata, locirepeats, nruns, burnin, record_interval)

