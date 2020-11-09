import subprocess
import sys
import csv


def main(argv):

    results = open('integration_test_results.csv', 'w')
    results_writter = csv.writer(results)

    for i in range(int(argv[0])):

        # run both scripts
        # EDIT FOR YOUR COMPUTER
        retR = subprocess.call(["C:/Program Files/R/R-4.0.2/bin/Rscript.exe", "--vanilla", "C:/Users/Heather/Dropbox/School/Senior/Junior Design/BayesianMicrosatellite/main.R"])
        retP = subprocess.call(["python3", "anti-malarial-MCMC-bayesian-algorithm/main.py"])

        #open output files
        p_prob_file = open('anti-malarial-MCMC-bayesian-algorithm/probability_of_recrudescence.csv', 'r')
        p_micro_file = open('anti-malarial-MCMC-bayesian-algorithm/microsatellite_correction.csv', 'r')
        r_prob_file = open('probability_of_recrudescence.csv', 'r')
        r_micro_file = open('microsatellite_correction.csv', 'r')
        #create results file

        p_prob_reader = csv.reader(p_prob_file)
        p_micro_reader = csv.reader(p_micro_file)
        r_prob_reader = csv.reader(r_prob_file)
        r_micro_reader = csv.reader(r_micro_file)

        # calculate stats for Probability of Recrudescence
        index = 0
        res = {'Num Samples':0, 'Total Diff': 0, 'Avg Diff':0, 'Min Diff':0, 'Max Diff':0}
        diffs = []
        lis = list(p_prob_reader)
        for row in r_prob_reader:
            if index == 0:
                results_row = ['Probability of Recrudescence Results: test run ' + str(i)]
            else:
                r_out = float(row[0])
                row2 = lis[index]
                p_out = float(row2[0])
                diff = abs(r_out-p_out)
                diffs.append(diff)
            index = index + 1
        results_writter.writerow(results_row)

        res['Total Diff'] = sum(diffs)
        res['Avg Diff'] = sum(diffs)/len(diffs)
        res['Min Diff'] = min(diffs)
        res['Max Diff'] = max(diffs)
        res['Num Samples'] = len(diffs)
        row1 = res.keys()
        results_writter.writerow(row1)
        row2 = res.values()
        results_writter.writerow(row2)
        results_writter.writerow([])


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
