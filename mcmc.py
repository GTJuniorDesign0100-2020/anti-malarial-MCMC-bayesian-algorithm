import numpy as np
import pandas as pd
import scipy.stats as sp_stats

from define_alleles import *
from calculate_frequencies import *
from recode_alleles import *
from switch_hidden import *
from findposteriorfrequencies import *

# TODO: Refactor this abomination. Then purge it with fire.
def onload(
    # global variables
    genotypedata_RR, additional_neutral, locirepeats, nruns, burnin, record_interval, jobname,
    ):
    # MOI = multiplicity of infection
    maxMOI = int(np.nanmax(  # Return array max, ignoring NaNs
        # NOTE: Assuming genotypedata_RR is a pandas dataframe
        # Split string like so: https://cmdlinetips.com/2018/06/how-to-split-a-column-or-column-names-in-pandas-and-get-part-of-it/
        # Gets the
        pd.to_numeric(genotypedata_RR.columns.str.split("_").str[1])
    ))

    # Get the unique Sample IDs in the dataset
    # NOTE: pd.unique used instead of np.unique to preserve ordering
    ids = pd.unique(
        genotypedata_RR[genotypedata_RR["Sample ID"].str.contains("Day 0")][
            "Sample ID"
        ].str.replace(" Day 0", "")
    )
    # Ditto, the unique loci for the set
    locinames = pd.unique(genotypedata_RR.columns[1:].str.split("_").str[0])

    nids = ids.size
    nloci = locinames.size

    maxalleles = 30

    k = np.repeat(maxalleles, nloci)

    alleles_definitions_RR = define_alleles(
        pd.concat([genotypedata_RR, additional_neutral]), locirepeats, k
    )

    ##### calculate MOI
    # TODO: What do MOI0 and MOIf stand for? Multiplicity of infection on first day/day of failure? What are they used for/mean?
    MOI0 = np.repeat(0, nids)
    MOIf = np.repeat(0, nids)
    for i, ID in enumerate(ids):
        for lociname in locinames:
            locicolumns = genotypedata_RR.columns.str.contains(f"{lociname}_")

            nalleles0 = np.count_nonzero(
                ~genotypedata_RR.loc[
                    genotypedata_RR["Sample ID"].str.contains(f"{ID} Day 0"), locicolumns
                ].isna()
            )
            nallelesf = np.count_nonzero(
                ~genotypedata_RR.loc[
                    genotypedata_RR["Sample ID"].str.contains(f"{ID} Day Failure"),
                    locicolumns,
                ].isna()
            )

            MOI0[i] = np.max([MOI0[i], nalleles0])
            MOIf[i] = np.max([MOIf[i], nallelesf])

    ##### define statevector (i.e. initial MCMC state)
    # TODO: What are each of these arrays used for? Should they not be initialized until closer to where they're actually used?
    alleles0 = np.zeros((nids, maxMOI * nloci))
    recoded0 = np.zeros((nids, maxMOI * nloci))
    hidden0 = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    recr0 = np.full_like(np.empty((nids, nloci)), np.nan)
    recr_repeats0 = np.full_like(
        np.empty((nids, nloci)), np.nan
    )  # number of times recrudescing allele is repeated on day 0
    recr_repeatsf = np.full_like(
        np.empty((nids, nloci)), np.nan
    )  # number of times recrudescing allele is repeated on day 0
    allelesf = np.zeros((nids, maxMOI * nloci))
    recodedf = np.zeros((nids, maxMOI * nloci))
    hiddenf = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    recrf = np.full_like(np.empty((nids, nloci)), np.nan)
    mindistance = np.zeros((nids, nloci))
    alldistance = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), np.nan)
    allrecrf = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), np.nan)
    classification = np.repeat(0, nids)

    ##### create state 0 (essentially, set NAN values to 0 and fill alleles with the appropriate initial/failure data)
    for i, locus in enumerate(locinames):
        # locicolumns code is duplicated from MOI calculations
        locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

        oldalleles = genotypedata_RR.loc[:, locicolumns].to_numpy()
        newalleles = np.copy(oldalleles)
        ncolumns = oldalleles.shape[1]
        for j in range(ncolumns):
            newalleles[:,j] = np.array(list(map(
                lambda x: recodeallele(alleles_definitions_RR[i].to_numpy(), oldalleles[x,j]),
                range(0, oldalleles.shape[0])
                )))

        # Set all nans in either array to 0
        oldalleles[np.isnan(oldalleles)] = 0
        oldalleles[np.isnan(newalleles)] = 0
        newalleles[np.isnan(newalleles)] = 0

        startColumn = maxMOI * i  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
        endColumnOldAllele = maxMOI * i + oldalleles.shape[1]
        endColumnNewAllele = maxMOI * i + newalleles.shape[1]
        alleles0[:, startColumn:endColumnOldAllele] = oldalleles[
            genotypedata_RR["Sample ID"].str.contains("Day 0"), :
        ]
        allelesf[:, startColumn:endColumnOldAllele] = oldalleles[
            genotypedata_RR["Sample ID"].str.contains("Day Failure"), :
        ]
        recoded0[:, startColumn:endColumnNewAllele] = newalleles[
            genotypedata_RR["Sample ID"].str.contains("Day 0"), :
        ]
        recodedf[:, startColumn:endColumnNewAllele] = newalleles[
            genotypedata_RR["Sample ID"].str.contains("Day Failure"), :
        ]

    ##### recode additional_neutral, if needed
    # TODO: What does recoding do? Why is it needed?
    # TODO: Seems to copy-paste much of the previous code section
    recoded_additional_neutral = None
    if additional_neutral.size > 0 and additional_neutral.shape[0] > 0:
        recoded_additional_neutral = np.zeros((additional_neutral.shape[0], maxMOI * nloci))
        for i, locus in enumerate(locinames):
            locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

            oldalleles = additional_neutral.loc[:, locicolumns].to_numpy()
            newalleles = np.copy(oldalleles)
            ncolumns = oldalleles.shape[1]

            for j in range(ncolumns):
                newalleles[:,j] = np.array(list(map(
                    lambda x: recodeallele(alleles_definitions_RR[i].to_numpy(), oldalleles[x,j]),
                    range(0, oldalleles.shape[0]))))
            newalleles[np.isnan(newalleles)] = 0
            oldalleles[np.isnan(oldalleles)] = 0

            oldalleles[newalleles == 0] = 0

            startColumn = maxMOI * (
                i - 1
            )  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
            endColumnOldAllele = maxMOI * (i - 1) + oldalleles.shape[1]
            recoded_additional_neutral[:, startColumn:endColumnOldAllele] = newalleles

    ## estimate frequencies
    frequencies_RR = calculate_frequencies3(
        pd.concat([genotypedata_RR, additional_neutral]), alleles_definitions_RR
    )

    ## assign random hidden alleles
    # TODO: Figure out what this code is overall trying to do (replace non-0 elements with random values? Why is it important that the values are assigned from frequencies_RR?)
    for i in range(nids):
        for j in range(nloci):
            # TODO: Code almost duplicated between top/bottom portions; refactor into single function? (needs 9 inputs: maxMOI, nids, nloci, MOIarray, alleles/recoded/hidden array, alleles_definitions_RR, frequencies_RR)
            # TODO: Start/end of what? The portion of the row w/ this locus information?
            start = maxMOI * j
            end = maxMOI * (j + 1)

            nalleles0 = np.count_nonzero(alleles0[i, start:end])
            nmissing0 = MOI0[i] - nalleles0

            # TODO: Rename "nonzero_indices" and "zero_indices"?
            whichnotmissing0 = np.arange(start, end)[
                np.where(alleles0[i, start : start + MOI0[i]] != 0)
            ]
            whichmissing0 = np.arange(start, end)[
                np.where(alleles0[i, start : start + MOI0[i]] == 0)
            ]

            # Sample to randomly initialize the alleles/hidden variables
            if nalleles0 > 0:
                hidden0[i, whichnotmissing0] = 0
            if nmissing0 > 0:
                newhiddenalleles0 = np.random.choice(
                    np.arange(
                        0, int(frequencies_RR[0][j])
                    ),  # Select from first row (count of how many probabilities they are)
                    size=nmissing0,
                    replace=True,
                    p=frequencies_RR[1][j, 0: int(frequencies_RR[0][j])]
                )  # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
                recoded0[i, whichmissing0] = newhiddenalleles0
                # calculate row means
                alleles0[i, whichmissing0] = np.mean(alleles_definitions_RR[j], axis=1)[
                    newhiddenalleles0
                ]  # hidden alleles get mean allele length
                hidden0[i, whichmissing0] = 1

            nallelesf = np.count_nonzero(allelesf[i, start:end])
            nmissingf = MOIf[i] - nallelesf

            # TODO: Rename "nonzero_indices" and "zero_indices"?
            whichnotmissingf = np.arange(start, end)[
                np.where(allelesf[i, start : start + MOIf[i]] != 0)
            ]
            whichmissingf = np.arange(start, end)[
                np.where(allelesf[i, start : start + MOIf[i]] == 0)
            ]

            if nallelesf > 0:
                hiddenf[i, whichnotmissingf] = 0
            if nmissingf > 0:
                newhiddenallelesf = np.random.choice(
                    np.arange(
                        0, int(frequencies_RR[0][j])
                    ),  # Select from first row (count of how many probabilities they are)
                    size=nmissingf,
                    replace=True,
                    p=frequencies_RR[1][j, 0 : int(frequencies_RR[0][j])]
                )  # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
                recodedf[i, whichmissingf] = newhiddenallelesf
                # calculate row means
                allelesf[i, whichmissingf] = np.mean(alleles_definitions_RR[j], axis=1)[
                    newhiddenallelesf
                ]  # hidden alleles get mean allele length
                hiddenf[i, whichmissingf] = 1

    ## initial estimate of q (probability of an allele being missed)
    # TODO: What is qq? Is there a better name for this?
    qq = np.nanmean(np.concatenate([hidden0, hiddenf]))

    ## initial estimate of dvect (likelihood of error in analysis)
    # TODO: What does dvect stand for?
    ranges = []
    for dataframe in alleles_definitions_RR:
        # Get the range (max-min) of the first "nloci" dataframes, then the max of all those
        ranges.append(dataframe.max().max() - dataframe.min().min())

    dvect = np.zeros(1 + int(round(max(ranges))))
    dvect[1] = 0.75
    dvect[2] = 0.2
    dvect[3] = 0.05

    ## randomly assign recrudescences/reinfections
    for i in range(nids):
        z = np.random.uniform(size=1)
        if z < 0.5:
            classification[i] = 1
        for j in range(
            nloci
        ):  # determine which alleles are recrudescing (for beginning, choose closest pair)
            allpossiblerecrud = np.stack(
                np.meshgrid(np.arange(MOI0[i]), np.arange(MOIf[i]))
            ).T.reshape(-1, 2)

            allele0_col_indices = maxMOI * j + allpossiblerecrud[:, 0]
            allelef_col_indices = maxMOI * j + allpossiblerecrud[:, 1]

            recrud_distances = np.abs(
                alleles0[i, allele0_col_indices] - allelesf[i, allelef_col_indices]
            )
            # rename to "closest_recrud_index"?
            closestrecrud = np.argmin(recrud_distances)

            mindistance[i, j] = recrud_distances[closestrecrud]
            alldistance[i, j, : recrud_distances.size] = recrud_distances

            allrecrf[i, j, : allpossiblerecrud.shape[0]] = recodedf[
                i, maxMOI * j + allpossiblerecrud[:, 1]
            ]
            recr0[i, j] = maxMOI * j + allpossiblerecrud[closestrecrud, 0]
            recrf[i, j] = maxMOI * j + allpossiblerecrud[closestrecrud, 1]

            recr_repeats0[i, j] = np.sum(
                recoded0[i, maxMOI * j : maxMOI * (j + 1)] == recoded0[i, int(recr0[i, j])]
            )  # TODO: int() only needed for stub
            recr_repeatsf[i, j] = np.sum(
                recodedf[i, maxMOI * j : maxMOI * (j + 1)] == recodedf[i, int(recrf[i, j])]
            )  # TODO: int() only needed for stub

    #### correction factor (reinfection)
    correction_distance_matrix = [] # for each locus, matrix of distances between each allele
    # TODO: Vectorize this (it seems fairly doable)
    for i in range(nloci):
        # Wrap mean call in "array" so we get a 2D array we can transpose (getting us a grid of distances, not just a 1D vector)
        distances = np.array([np.mean(alleles_definitions_RR[i], axis=1)])
        distance_combinations = np.abs(distances.T - distances)
        correction_distance_matrix.append(distance_combinations)

    num_saved_records = int((nruns - burnin) / record_interval)

    state_classification = np.full_like(np.empty((nids, num_saved_records)), np.nan)
    state_alleles0 = np.full_like(
        np.empty((nids, maxMOI * nloci, num_saved_records)), np.nan
    )
    state_allelesf = np.full_like(
        np.empty((nids, maxMOI * nloci, num_saved_records)), np.nan
    )
    state_parameters = np.full_like(np.empty((2 + 2 * nloci, num_saved_records)), np.nan)

    dposterior = 0.75


    def runmcmc(iteration, dvect, classification):
        # propose new classification
        likelihoodratio = np.zeros(nids)
        # TODO: Finish vectorizing this
        for x in range(nids):
            # id mean for what?
            id_means = np.zeros(nloci)
            for y in range(nloci):
                id_means[y] = np.nanmean(
                    dvect[np.round(alldistance[x, y, :][~np.isnan(alldistance[x, y, : ])]).astype(int)]
                    # Should get an array of maxMOI**2 sums
                    / np.sum(
                        # TODO: Make sure multiplications are down the right axis (I believe each element in the frequencies_RR 1D vector should multiply across 1 dvect row)
                        # Double-transpose to multiply across rows, not columns
                        (frequencies_RR[1][y, :int(frequencies_RR[0][y])]
                        * dvect[
                            correction_distance_matrix[y][
                                :,
                                allrecrf[x, y, :maxMOI**2][~np.isnan(allrecrf[x, y, :maxMOI**2])].astype(int),
                            ].astype(int)
                        ].T).T,
                        axis=0, # TODO: Verify it's the right axis?
                    ),
                )
            likelihoodratio[x] = np.exp(np.sum(np.log(id_means)))

        z = np.random.uniform(size=nids)
        newclassification = classification
        newclassification[np.logical_and(classification == 0, z < likelihoodratio)] = 1
        newclassification[np.logical_and(classification == 1, z < 1 / likelihoodratio)] = 0
        classification = newclassification

        # propose new hidden states
        # TODO: What does switch_hidden do? Is it entirely side effects? (Also,: can't run this yet, still waiting on implementation)
        for i in range(nids):
            switch_hidden(i, hidden0, hiddenf)

        # propose q (beta distribution is conjugate distribution for binomial process)
        q_prior_alpha = 0
        q_prior_beta = 0
        q_posterior_alpha = (
            q_prior_alpha + np.nansum(hidden0 == 1) + np.nansum(hiddenf == 1)
        )
        q_posterior_beta = q_prior_beta + np.nansum(hidden0 == 0) + np.nansum(hiddenf == 0)
        if q_posterior_alpha == 0:
            q_posterior_alpha = 1
        if q_posterior_beta == 0:  # TODO: Added this due to numpy warning, possibly remove?
            q_posterior_beta = 1
        qq = np.random.beta(q_posterior_alpha, q_posterior_beta)

        #  update dvect (approximate using geometric distribution)
        # only if there is at least 1 recrudescent infection
        if np.sum(classification == 1) >= 1:
            d_prior_alpha = 0
            d_prior_beta = 0
            d_posterior_alpha = d_prior_alpha + mindistance[classification == 1, :].size
            d_posterior_beta = d_prior_beta + np.sum(
                np.round(mindistance[classification == 1, :])
            )
            if d_posterior_beta == 0:
                d_posterior_beta = np.sum(mindistance[classification == 1, :])
            if (
                d_posterior_beta == 0
            ):  ## algorithm will get stuck if dposterior is allowed to go to 1 (TODO: Wait, so why is it setting d_posterior_beta to 1??)
                d_posterior_beta = 1

            dposterior = np.random.beta(d_posterior_alpha, d_posterior_beta)
            dvect = dposterior * (np.array(1 - dposterior) ** np.arange(0, dvect.size))
            dvect = dvect / np.sum(dvect)

        # update frequencies
        # remove recrudescing alleles from calculations
        tempdata = np.copy(recoded0)
        recrudescent_alleles = np.where(classification == 1)
        tempdata[recrudescent_alleles, recr0[recrudescent_alleles, :]] = 0
        tempdata = pd.concat([tempdata, recodedf])
        for i in range(nloci):
            findposteriorfrequencies(x, pd.concat([tempdata, recoded_additional_neutral]))

        # record state
        if iteration > burnin and iteration % record_interval == 0:
            print(iteration)
            record_index = (iteration - burnin) / record_interval
            state_classification[:, record_index] = classification
            state_alleles0[:, :, record_index] = alleles0
            state_allelesf[:, :, record_index] = allelesf
            state_parameters[0, record_index] = qq
            state_parameters[1, record_index] = dposterior
            state_parameters[2 : (2 + nloci), record_index] = frequencies_RR[1].max(axis=0)
            state_parameters[2 + nloci : (2 + 2 * nloci), record_index] = np.sum(
                frequencies_RR[1, :nloci, :] ** 2
            )


    for i in range(nruns):
        runmcmc(i, dvect, classification)

    ## make sure no NAs in result matrices
    state_parameters = state_parameters[:, ~np.isnan(np.sum(state_parameters, axis=1))]
    state_classification = state_classification[
        :, ~np.isnan(np.sum(state_classification, axis=1))
    ]

    modealleles = np.zeros((2 * nids, maxMOI * nloci))
    for i in range(nids):
        for j in range(nloci):
            modealleles[2 * i, j * maxMOI : (j + 1) * maxMOI] = sp_stats.mode(
                state_alleles0[i, j * maxMOI : (j + 1) * maxMOI, :], axis=1
            )[0].ravel()

            modealleles[2 * i + 1, j * maxMOI : (j + 1) * maxMOI] = sp_stats.mode(
                state_allelesf[i, j * maxMOI : (j + 1) * maxMOI, :], axis=1
            )[0].ravel()

    # TODO: Combined what?
    temp_combined = np.repeat(np.mean(state_classification, axis=0)[:nids], 2)
    outputmatrix = np.concatenate((temp_combined, modealleles), axis=1)
    outputmatrix_columns = [
        [f"{locus}_{i+1}" for i in range(maxMOI)] for locus in locinames
    ]
    outputmatrix_columns.insert(0, "Prob Rec")
    outputmatrix = pd.DataFrame(outputmatrix, columns=outputmatrix_columns)
    outputmatrix.to_csv(f"{jobname}_posterior.csv")

    # summary statistics of parameters
    pd.DataFrame(state_parameters).to_csv(f"{jobname}_state_parameters.csv")

    summary_statisticsmatrix = np.concatenate(
        (
            np.mean(state_parameters, axis=1).reshape(-1, 1),
            np.quantile(state_parameters, (0.25, 0.75), axis=1).T,
        ),
        axis=1,
    )
    summary_statisticsmatrix = np.concatenate(
        (
            summary_statisticsmatrix,
            np.append(
                np.quantile(state_parameters[2 + nloci :, :], (0.25, 0.75)),
                np.mean(state_parameters[2 + nloci :, :]),
            ).reshape(1, -1),
        )
    )
    summary_statisticsmatrix = np.array(
        [
            f"{summary_statisticsmatrix[i,0]:.2f} ({summary_statisticsmatrix[i,1]:.2f}, {summary_statisticsmatrix[i,2]:.2f})"
            for i in range(summary_statisticsmatrix.shape[0])
        ]
    )
    summary_statisticsmatrix_df = pd.DataFrame(
        summary_statisticsmatrix,
        index=["q", "d", *locinames.tolist(), *locinames.tolist(), "Mean diversity"],
    )
    summary_statisticsmatrix_df.to_csv(f"{jobname}_summarystatistics.csv")

    return state_classification, state_parameters, ids
