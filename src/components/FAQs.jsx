import React from 'react';



export default class FAQs extends React.Component {
  constructor(props) {
      super(props);
  }

  render() {
      return (
        <div class='FAQData'>
            <h1>FAQs</h1>
                <h2>How do I format the input data</h2>
                    <p>The input data should be an Excel file with 2 worksheets. The first worksheet holds the paired Day 0 – Day of Failure data, while the second worksheet holds additional background (unpaired) data. Additional background samples are not necessary for the algorithm to run, but the precision of the estimates will be improved the more background genotypes are included. For both datasets, rows are samples, and columns correspond to microsatellite loci. The columns are indexed by the number of alleles present for each locus. For example, the column “Polya_2” contains the 2nd allele for the “Polya” locus. The values of the cells contain the measured sample lengths for each sample. Paired samples should be labeled as SAMPLEID_D0 and SAMPLEID_DX, where SAMPLEID is the unique sample identification number, and X is the Day of Failure (eg. 28). The “Site” column should contain the name of the sentinel site where samples were collected; all samples from the same site are analyzed together by the algorithm. <br/><br/>
                    The “loci repeats” input parameter represents the vector of sizes of the tandem repeats for each marker. This is an inherent property of each microsatellite marker and is known a priori. For example, the PolyA microsatellite marker consists of repeated trimers of the form ATT, and its loci repeat value is therefore 3.<br/><br/>
                    A sample of the input data format can be found here: <a href="Angola_2017_example.xlsx" download>Angola Sample Data</a></p>
                <h2>What does the output mean?</h2>
                    <p>This application outputs Comma Separated Values (.csv) files that can be opened in Microsoft Excel or Google Sheets. The file probability_of_recrudescence.csv contains the posterior probability of recrudescence for each episode of recurrent parasitemia. This value represents the weight of the statistical evidence that the genotype data indicate that at least one strain on Day of Failure was also present on Day 0. It ranges from 0 (no evidence of recrudescence) to 1 (full evidence of recrudescence). Indeterminate values indicate that there is some evidence, but not full certainty, for recrudescence. For example, a posterior probability of recrudescence of 0.75 can be interpreted as a 75% likelihood that the episode of recurrent parasitemia represents a recrudescence and a 25% chance that it is a reinfection.
                    <br/><br/>Further statistics and information regarding the parameters of the tests are available in a ZIP archive. These additional files contain estimates of the hidden, unobserved states (files ending in “_posterior.csv”), and summary statistics (files ending in “_summarystatistics.csv”) for each site. Researchers can analyze these results to gain further insight into the properties inferred for each site. For example, estimates of diversity can allow researchers to identify which markers are the most diverse and thus the most informative.</p>
                <h2>How many iterations should I run?</h2>
                    <p>The minimum number of iterations will depend on the size of the sample set. Larger sample sets will need more iterations for the algorithm to converge. A good starting point is 10,000-100,000 iterations. The simplest way to judge convergence is to run the algorithm a few times on the same dataset. If the outputs are consistent from run to run, then the algorithm has converged. Strategies for rigorously assessing convergence for this kind of algorithm are still a subject of <a href="https://www.tandfonline.com/doi/abs/10.1080/01621459.1996.10476956">ongoing research in the statistical community</a>.</p>
                <h2>How do I use the output to calculate my corrected efficacy?</h2>
                    <p>There are several options to use the posterior probability of recrudescence to estimate the corrected efficacy (in decreasing order of preference):
                        <br/>1. Generate a sample of recrudescence/reinfection calls from the posterior probability of recrudescence and use the Kaplan-Meier estimator to generate an empiric sample of the corrected efficacy. A sample R script using this approach can be found here.
                        <br/>2. Use the sum of the posterior probability of recrudescence as the estimate of the true number of recrudescences. For example, a study with 70 treatment successes, 30 episodes of recurrent parasitemia, and a sum of the posterior probability of recrudescence of 12.5, means that the best estimate is that there were 12.5 recrudescences (and 17.5 reinfections) amongst the 30 episodes of recurrent parasitemia. The corrected efficacy would then be 70 / (70 + 12.5) or 84.8%.
                        <br/>3. Dichotomize the episodes of recurrent parasitemia as recrudescences or reinfections using a threshold. Research shows that the most unbiased threshold for this method would be 0.5, i.e. any patient with posterior probability of recrudescence > 0.5 would be considered a recrudescence and any patient with posterior probability of recrudescence ≤ 0.5 would be considered a reinfection.</p>
                <h2>Has this tool been validated?</h2>
                    <p>The algorithm was first developed in 2015, and the scientific report describing it can be found <a href="https://aac.asm.org/content/59/10/6096">here</a>. Since then, the algorithm has been <a href="https://aac.asm.org/content/64/4/e01517-19">validated by modelling studies</a> and is recommended by the Centers for Disease Control and Prevention as the standard of practice for analysis of microsatellite datasets.</p>
                <h2>What are Bayesian statistics</h2>
                    <p>Bayesian statistical methods represent one of the primary schools of thought in modern statistics. They are particularly well suited for complex statistical problems with multiple parameters and unobserved states. They have a <a href="https://www.annualreviews.org/doi/abs/10.1146/annurev.pu.16.050195.000323">long history of use in biology and epidemiology</a>.</p>
                <h2>Why are the results different every time I run the algorithm?</h2>
                    <p>The statistical problem of inferring recrudescence and reinfection using microsatellite data with unobserved states is too complex for standard statistical methods. For computationally intractable problems like this, stochastic algorithms are necessary. Stochastic algorithms approximate parameter estimates, but results will vary slightly from run to run by definition. An accessible introduction to Monte Carlo Markov Chain algorithms, the type of algorithm implemented here, <a href="https://academic.oup.com/ije/article/42/2/627/738896">can be found here</a>.</p>
                <h2>Can I run this locally?</h2>
                    <p>Yes! If you want to run this locally or develop additions, changes, or improvements, you can find all the code for the application on GitHub, <a href="https://github.com/GTJuniorDesign0100-2020/anti-malarial-MCMC-bayesian-algorithm">here</a>! The application runs a ReactJS front end, a flask API, and a Python back end, and can be configured locally to run without the front end if desired. To do so, update input variables and run the main.py file. </p>
        </div>
    );
  }
}
