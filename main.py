import os
import numpy as np

# set working directory 
# os.chdir("\\\\cdc.gov\\private\\M326\\wif7\\Microsatellite\\GATech\\MRE");

#import gtools;


# user inputs
inputfile = "Angola2017_example.xlsx";
locirepeats = np.array([[2,2,3,3,3,3,3]]);
nruns = 1000;

# call script to import data
os.system("python Import_Microsatellite_Data.py");


# calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
record_interval = np.ceil(nruns / 1000);
burnin = np.ceil(nruns * 0.25);
os.system("python run_all_arms.py");

