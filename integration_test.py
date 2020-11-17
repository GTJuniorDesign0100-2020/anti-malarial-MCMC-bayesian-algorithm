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
        #retR = subprocess.call(["C:/Program Files/R/R-4.0.3/bin/Rscript.exe", "--vanilla", "C:/Users/Heather/JD/BayesianMicrosatellite/main.r"])
        retP = subprocess.call(["python3", "anti-malarial-MCMC-bayesian-algorithm/main.py"])

        por_data_r.append(get_por('r'))
        por_data_p.append(get_por('p'))

        mc_data_r.append(get_mc('r'))
        mc_data_p.append(get_mc('p'))


    sum_p = np.array(por_data_p).T.tolist()
    p_conv = []
    for x in range(len(sum_p[1])):
        if x != 0:
            fl = [float(el) for el in sum_p[1][x]]
            p_conv.append(sum(fl)/(len(sum_p[1])-1))

    sum_r = np.array(por_data_r).T.tolist()
    r_conv = []
    fl = []
    for x in range(len(sum_r[1])):
        if x != 0:
            fl = [float(el) for el in sum_r[1][x]]
            r_conv.append(sum(fl)/(len(sum_r[1])-1))

#    print(p_conv)
#    print(r_conv)

    diff = []
    for j in range(len(p_conv)):
        diff.append(abs(p_conv[j]-r_conv[j]))

    results_row = ["Probability of Recrudescence Results: After " + str(i+1) + " iterations (of 1000 runs):"]
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

    res = {'Num Samples':0, 'Num Itterations': 0, 'Total Compairisons':0, 'Similarities':0, 'Differences':0, 'Percent Correct':0}
    p = np.array(mc_data_p).T.tolist()
    r = np.array(mc_data_r).T.tolist()
    res['Num Itterations'] = len(r)
    res['Num Samples'] = len(p[0])
    num_samp = res['Num Samples'] - 1
    for index in range(len(p)):
        if index != 0 and index != len(p)-1:
            r_values = r[index + 1][1:]
            p_values = p[index + 1][1:]
            r_values = [float(x[0]) for x in r_values]
            p_values = [float(x[0]) for x in p_values]
            for k in range(len(r_values)):
                res['Total Compairisons'] += 1
                if r_values[k] == p_values[k]:
                    res['Similarities'] += 1
                else:
                    res['Differences'] += 1
    results_row = ["Microsatellite Correction Results: After " + str(i+1) + " iterations (of 1000 runs):"]
    res['Percent Correct'] = res['Similarities']/res['Total Compairisons']
    row1 = res.keys()
    results_writter.writerow(results_row)
    results_writter.writerow(row1)
    row2 = res.values()
    results_writter.writerow(row2)




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



if __name__=="__main__":
    main(sys.argv[1:])
