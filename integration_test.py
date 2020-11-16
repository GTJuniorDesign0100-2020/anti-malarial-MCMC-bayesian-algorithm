import subprocess
import sys
import csv
import numpy as np


def main(argv):

    results = open('integration_test_results.csv', 'w')
    results_writter = csv.writer(results)

    por_data_r = []
    mc_data_r = []

    por_data_p = []
    mc_data_p = []

    ittr = int(argv[0])

    for i in range(ittr):

        # run both scripts
        # EDIT FOR YOUR COMPUTER
        retR = subprocess.call(["C:/Program Files/R/R-4.0.3/bin/Rscript.exe", "--vanilla", "C:/Users/Heather/JD/BayesianMicrosatellite/main.r"])
        retP = subprocess.call(["python3", "anti-malarial-MCMC-bayesian-algorithm/main.py"])

        por_data_r.append(get_por('r'))
        por_data_p.append(get_por('p'))


    sum_p = np.array(por_data_p).T.tolist()
    p_conv = []
    for x in range(len(sum_p[1])):
        if x != 0:
            fl = [float(el) for el in sum_p[1][x]]
            p_conv.append(sum(fl)/(len(sum_p[1])-1))

    sum_r = np.array(por_data_r).T.tolist()
    r_conv = []
    for x in sum_r[0]:
        fl = [float(el) for el in x]
        r_conv.append(sum(fl)/len(x))

    diff = []
    for j in range(len(p_conv)):
        diff.append(abs(p_conv[j]-r_conv[j]))

    results_row = ["After " + str(i+1) + " iterations (of 1000 runs):"]
    results_writter.writerow(results_row)
    results_row = ["R code Converged Too:"] + r_conv
    results_writter.writerow(results_row)
    results_row = ["Python code Converged Too:"] + p_conv
    results_writter.writerow(results_row)
    results_row = ["The Total differences are:"] + diff
    results_writter.writerow(results_row)

    res = {'Num Samples':0, 'Total Diff': 0, 'Avg Diff':0, 'Min Diff':0, 'Max Diff':0}
    res['Total Diff'] = sum(diff)
    res['Avg Diff'] = sum(diff)/len(diff)
    res['Min Diff'] = min(diff)
    res['Max Diff'] = max(diff)
    res['Num Samples'] = len(diff)

    row1 = res.keys()
    results_writter.writerow(row1)
    row2 = res.values()
    results_writter.writerow(row2)
    results_writter.writerow([])




def get_por(t):
    #open output files
    if t == 'p':
        prob_file = open('anti-malarial-MCMC-bayesian-algorithm/probability_of_recrudescence.csv', 'r')

    if t == 'r':
        prob_file = open('probability_of_recrudescence.csv', 'r')

    prob_reader = csv.reader(prob_file)
    return list(prob_reader)


def get_mc(t):
    #open output files
    if t == 'r':
        micro_file = open('microsatellite_correction.csv', 'r')

    if t == 'p':
        micro_file = open('anti-malarial-MCMC-bayesian-algorithm/microsatellite_correction.csv', 'r')

    micro_reader = csv.reader(micro_file)
    return list(micro_reader)


    # calculate stats for microsatellite correction
    """
    index = 0
    res = {'Num Samples':0, 'Num Itterations': 0, 'Total Compairisons':0, 'Similarities':0, 'Differences':0, 'Percent Correct':0}
    results_row = []
    lis = list(p_micro_reader)
    for row in r_micro_reader:
        if index == 0:
            results_row = ['Microsatellite Correction Results: test run ' + str(i), 'Num Similar', 'Num Differnt']
        else:
            num_ittr = len(row) - 1
            res['Num Itterations'] = num_ittr
            rowp = lis[index + 39]
            correct = 0
            incorrect = 0
            for ittr in range(num_ittr):
                r_out = bool(row[ittr + 1])
                p_out = bool(rowp[ittr])
                if r_out == p_out:
                    correct += 1
                    res['Similarities'] += 1
                else:
                    incorrect += 1
                    res['Differences'] += 1
            results_row = [str(lis[index][0]), str(correct), str(incorrect)]
        index = index + 1
        results_writter.writerow(results_row)
    res = {'Num Samples':0, 'Num Itterations': 0, 'Total Compairisons':0, 'Similarities':0, 'Differences':0, 'Percent Correct':0}
    res['Num Samples'] = index
    res['Total Compairisons'] = index*res['Num Itterations']
    res['Percent Correct'] = res['Similarities']/res['Total Compairisons']
    row1 = res.keys()
    results_writter.writerow(row1)
    row2 = res.values()
    results_writter.writerow(row2)
    """


if __name__=="__main__":
    main(sys.argv[1:])
