import subprocess
import sys
import csv
import numpy as np

'''
This integration test algorithm tests the python algorithm output against the
original R script outputs. Specifically, it checks microsatellite_correction.csv
and probability_of_recrudescence.csv. It takes 1 argument that is a simple
integer and represents the number of times you run the algorithm (each doing
1000 iterations), so for example running with 'python integration_test.py 5'
will run main.py 5 times. At the end, it outputs a csv file with the results
in a file called integration_test_results.csv. To run this test, update
lines that have 'UPDATE HERE' commented before.

'''


def main(argv):

    results = open('integration_test_results.csv', 'w')
    results_writer = csv.writer(results)

    por_data_r = []
    mc_data_r = []

    por_data_p = []
    mc_data_p = []

    num_runs = int(argv[0])

    for i in range(num_runs):

        # run both scripts
        # UPDATE HERE : point to your local R executable file and Baysian algorithm in R
        retR = subprocess.call(["C:/Program Files/R/R-4.0.3/bin/Rscript.exe", "--vanilla", "C:/Users/Heather/JD/BayesianMicrosatellite/main.r"])
        # UPDATE HERE : point to your local python algorithm
        retP = subprocess.call(["python3", "anti-malarial-MCMC-bayesian-algorithm/main.py"])

        por_data_r.append(get_probability_of_recrudescence_file('r'))
        por_data_p.append(get_probability_of_recrudescence_file('p'))

        mc_data_r.append(get_microsatellite_correction_file('r'))
        mc_data_p.append(get_microsatellite_correction_file('p'))


    sum_p = np.array(por_data_p).T.tolist()
    p_conv = []
    #sum_p[1] represents the data from the python files
    for row in range(len(sum_p[1])):
        if row != 0:
            #since row is a list of strings we must convert to floats
            row_float = [float(el) for el in sum_p[1][row]]
            p_conv.append(sum(row_float)/(len(sum_p[1])-1))

    sum_r = np.array(por_data_r).T.tolist()
    r_conv = []
    fl = []
    #sum_r[1] represents the data from the R files
    for row in range(len(sum_r[1])):
        if row != 0:
            row_float = [float(el) for el in sum_r[1][row]]
            r_conv.append(sum(row_float)/(len(sum_r[1])-1))

    diff = []
    for j in range(len(p_conv)):
        diff.append(abs(p_conv[j]-r_conv[j]))

    results_row = ["Probability of Recrudescence Results: After " + str(i+1) + " iterations (of 1000 runs):"]
    results_writer.writerow(results_row)
    results_row = ["R code Converged Too:"] + r_conv
    results_writer.writerow(results_row)
    results_row = ["Python code Converged Too:"] + p_conv
    results_writer.writerow(results_row)
    results_row = ["The Total differences are:"] + diff
    results_writer.writerow(results_row)

    res = {'Num Samples':0, 'Total Diff': 0, 'Avg Diff':0, 'Min Diff':0, 'Max Diff':0}
    res['Total Diff'] = sum(diff)
    res['Avg Diff'] = sum(diff)/len(diff)
    res['Min Diff'] = min(diff)
    res['Max Diff'] = max(diff)
    res['Num Samples'] = len(diff)

    row1 = res.keys()
    results_writer.writerow(row1)
    row2 = res.values()
    results_writer.writerow(row2)
    results_writer.writerow([])

    res = {'Num Samples':0, 'Num Iterations': 0, 'Total Comparisons':0, 'Similarities':0, 'Differences':0, 'Percent Correct':0}
    p = np.array(mc_data_p).T.tolist()
    r = np.array(mc_data_r).T.tolist()
    res['Num Iterations'] = len(r)
    res['Num Samples'] = len(p[0])
    num_samp = res['Num Samples'] - 1
    for index in range(len(p)):
        if index != 0 and index != len(p)-1:
            r_values = r[index + 1][1:]
            p_values = p[index + 1][1:]
            r_values = [float(x[0]) for x in r_values]
            p_values = [float(x[0]) for x in p_values]
            for k in range(len(r_values)):
                res['Total Comparisons'] += 1
                if r_values[k] == p_values[k]:
                    res['Similarities'] += 1
                else:
                    res['Differences'] += 1
    results_row = ["Microsatellite Correction Results: After " + str(i+1) + " iterations (of 1000 runs):"]
    res['Percent Correct'] = res['Similarities']/res['Total Comparisons']
    row1 = res.keys()
    results_writer.writerow(results_row)
    results_writer.writerow(row1)
    row2 = res.values()
    results_writer.writerow(row2)




def get_probability_of_recrudescence_file(t):
    #open output files
    if t == 'p':
        # UPDATE HERE: to point to the probability_of_recrudescence.csv outputted by the python script
        prob_file = open('anti-malarial-MCMC-bayesian-algorithm/probability_of_recrudescence.csv', 'r')

    if t == 'r':
        # UPDATE HERE: to point to the probability_of_recrudescence.csv outputted by the R script
        prob_file = open('probability_of_recrudescence.csv', 'r')

    prob_reader = csv.reader(prob_file)
    return list(prob_reader)


def get_microsatellite_correction_file(t):
    #open output files
    if t == 'r':
        # UPDATE HERE: to point to the microsatellite_correction.csv outputted by the python script
        micro_file = open('microsatellite_correction.csv', 'r')

    if t == 'p':
        # UPDATE HERE: to point to the microsatellite_correction.csv outputted by the R script
        micro_file = open('anti-malarial-MCMC-bayesian-algorithm/microsatellite_correction.csv', 'r')

    micro_reader = csv.reader(micro_file)
    return list(micro_reader)



if __name__=="__main__":
    main(sys.argv[1:])
