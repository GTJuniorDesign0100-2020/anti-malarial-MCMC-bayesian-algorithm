"""
A series of tests to make sure mcmc.py's functionality is equivalent to mcmc.r's

TODO: Use an actual testing framework?
"""

import numpy as np
import pandas as pd

"""
Full example genotype_RR dataframe in R:

               Sample.ID X313_1 X313_2 X313_3 X383_1 X383_2 X383_3 TA1_1 TA1_2
1        BQ17-269_ Day 0  223.4     NA     NA  103.7  140.3     NA 169.0    NA
2  BQ17-269_ Day Failure  225.5     NA     NA  124.1  162.3     NA 178.0    NA
3        BD17-040_ Day 0  229.3     NA     NA  139.1     NA     NA 175.0    NA
4  BD17-040_ Day Failure  243.6     NA     NA  139.0     NA     NA 175.0    NA
5        BD17-083_ Day 0  221.3     NA     NA  151.6     NA     NA 165.2 180.6
6  BD17-083_ Day Failure  231.5     NA     NA  145.0     NA     NA 162.7 165.6
7        BD17-085_ Day 0  261.7     NA     NA  122.6     NA     NA 162.7    NA
8  BD17-085_ Day Failure  239.6     NA     NA  124.0     NA     NA 162.7    NA
9        BD17-087_ Day 0  261.7     NA     NA  124.0   89.0     NA 162.7    NA
10 BD17-087_ Day Failure  245.3     NA     NA  123.6  140.2     NA 163.0 184.0
11       BD17-090_ Day 0  238.4     NA     NA   87.0  123.9     NA 180.8    NA
12 BD17-090_ Day Failure  219.4  249.7     NA   87.0  123.6     NA 165.8 177.9
   TA1_3 POLYA_1 POLYA_2 POLYA_3 POLYA_4 POLYA_5 PFPK2_1 PFPK2_2 PFPK2_3
1     NA   167.3      NA      NA      NA      NA   149.8   158.0   171.2
2     NA   171.0      NA      NA      NA      NA   162.0   159.0      NA
3     NA   144.0      NA      NA      NA      NA   171.0   168.0      NA
4     NA   164.0      NA      NA      NA      NA   171.0      NA      NA
5     NA   170.2      NA      NA      NA      NA   168.2   177.3   180.3
6     NA   104.8   112.3      NA      NA      NA   168.0   162.0      NA
7     NA   167.1      NA      NA      NA      NA   194.7      NA      NA
8     NA   173.0      NA      NA      NA      NA   189.0      NA      NA
9     NA   167.3      NA      NA      NA      NA   186.0   191.9   194.8
10    NA   155.0      NA      NA      NA      NA   168.2   174.2   177.2
11    NA   151.3      NA      NA      NA      NA   161.9      NA      NA
12    NA   151.3   164.0      NA      NA      NA   161.9   168.1      NA
   PFPK2_4 X2490_1 X2490_2 TA109_1 TA109_2 TA109_3 TA109_4
1       NA    81.6      NA   164.6   175.0      NA      NA
2       NA    81.8      NA   161.6   172.5      NA      NA
3       NA    82.0      78   149.6   159.9      NA      NA
4       NA    82.0      NA   175.7   187.2      NA      NA
5       NA    81.8      NA   164.2   175.0      NA      NA
6       NA    79.0      NA   149.4   160.0      NA      NA
7       NA    78.6      NA   147.9   158.8      NA      NA
8       NA    88.0      NA   170.0   181.4      NA      NA
9       NA    78.6      NA   148.6   159.9      NA      NA
10      NA    82.0      NA   152.4   163.1      NA      NA
11      NA    81.7      NA   148.8   160.4      NA      NA
12      NA    81.7      NA   163.2   175.3      NA      NA
"""

# Hardcode this in for now
genotypedata_RR = pd.DataFrame(
    {
        "Sample.ID": [
            "BQ17-269_ Day 0",
            "BQ17-269_ Day Failure",
            "BD17-040_ Day 0",
            "BD17-040_ Day Failure",
            "BD17-083_ Day 0",
            "BD17-083_ Day Failure",
            "BD17-085_ Day 0",
            "BD17-085_ Day Failure",
            "BD17-087_ Day 0",
            "BD17-087_ Day Failure",
            "BD17-090_ Day 0",
            "BD17-090_ Day Failure",
        ],
        "X313_1": [
            223.4,
            225.5,
            229.3,
            243.6,
            221.3,
            231.5,
            261.7,
            239.6,
            261.7,
            245.3,
            238.4,
            219.4,
        ],
        "X313_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            249.7,
        ],
        "X313_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "X383_1": [
            103.7,
            124.1,
            139.1,
            139.0,
            151.6,
            145.0,
            122.6,
            124.0,
            124.0,
            123.6,
            87.0,
            87.0,
        ],
        "X383_2": [
            140.3,
            162.3,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            89.0,
            140.2,
            123.9,
            123.6,
        ],
        "X383_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA1_1": [
            169.0,
            178.0,
            175.0,
            175.0,
            165.2,
            162.7,
            162.7,
            162.7,
            162.7,
            163.0,
            180.8,
            165.8,
        ],
        "TA1_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            180.6,
            165.6,
            np.nan,
            np.nan,
            np.nan,
            184.0,
            np.nan,
            177.9,
        ],
        "TA1_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_1": [
            167.3,
            171.0,
            144.0,
            164.0,
            170.2,
            104.8,
            167.1,
            173.0,
            167.3,
            155.0,
            151.3,
            151.3,
        ],
        "POLYA_2": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            112.3,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            164.0,
        ],
        "POLYA_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "POLYA_5": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "PFPK2_1": [
            149.8,
            162.0,
            171.0,
            171.0,
            168.2,
            168.0,
            194.7,
            189.0,
            186.0,
            168.2,
            161.9,
            161.9,
        ],
        "PFPK2_2": [
            158.0,
            159.0,
            168.0,
            np.nan,
            177.3,
            162.0,
            np.nan,
            np.nan,
            191.9,
            174.2,
            np.nan,
            168.1,
        ],
        "PFPK2_3": [
            171.2,
            np.nan,
            np.nan,
            np.nan,
            180.3,
            np.nan,
            np.nan,
            np.nan,
            194.8,
            177.2,
            np.nan,
            np.nan,
        ],
        "PFPK2_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "X2490_1": [
            81.6,
            81.8,
            82.0,
            82.0,
            81.8,
            79.0,
            78.6,
            88.0,
            78.6,
            82.0,
            81.7,
            81.7,
        ],
        "X2490_2": [
            np.nan,
            np.nan,
            78,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA109_1": [
            164.6,
            161.6,
            149.6,
            175.7,
            164.2,
            149.4,
            147.9,
            170.0,
            148.6,
            152.4,
            148.8,
            163.2,
        ],
        "TA109_2": [
            175.0,
            172.5,
            159.9,
            187.2,
            175.0,
            160.0,
            158.8,
            181.4,
            159.9,
            163.1,
            160.4,
            175.3,
        ],
        "TA109_3": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        "TA109_4": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
    }
)


def test_max_MOI():
    maxMOI = np.nanmax(  # Return array max, ignoring NaNs
        # NOTE: Assuming genotypedata_RR is a pandas dataframe
        # Split string like so: https://cmdlinetips.com/2018/06/how-to-split-a-column-or-column-names-in-pandas-and-get-part-of-it/
        # Gets the
        pd.to_numeric(genotypedata_RR.columns.str.split("_").str[1])
    )
    expected_MOI = 5
    assert maxMOI == expected_MOI, f"{maxMOI} (expected {expected_MOI})"


def test_getting_ids():
    expected = np.unique(
        ["BQ17-269_", "BD17-040_", "BD17-083_", "BD17-085_", "BD17-087_", "BD17-090_"]
    )

    ids = np.unique(
        genotypedata_RR[genotypedata_RR["Sample.ID"].str.contains("Day 0")][
            "Sample.ID"
        ].str.replace(" Day 0", "")
    )

    assert np.array_equal(ids, expected), f"{ids} (expected {expected})"


def test_getting_locinames():
    expected = np.unique(["X313", "X383", "TA1", "POLYA", "PFPK2", "X2490", "TA109"])

    locinames = np.unique(genotypedata_RR.columns[1:].str.split("_").str[0])

    assert np.array_equal(locinames, expected), f"{locinames} (expected {expected})"


def test_calculate_MOI():
    # NOTE: These are different orderings than the original (possibly np.unique changes ordering?); I THINK this is okay, but should ask Mat
    expected_MOI0 = np.array([2, 3, 2, 3, 2, 3])
    expected_MOIf = np.array([2, 2, 2, 3, 2, 2])

    ids = np.unique(
        ["BQ17-269_", "BD17-040_", "BD17-083_", "BD17-085_", "BD17-087_", "BD17-090_"]
    )
    locinames = np.unique(["X313", "X383", "TA1", "POLYA", "PFPK2", "X2490", "TA109"])

    nids = ids.size
    nloci = locinames.size

    MOI0 = np.repeat(0, nids)
    MOIf = np.repeat(0, nids)
    for i, ID in enumerate(ids):
        for lociname in locinames:
            locicolumns = genotypedata_RR.columns.str.contains(f"{lociname}_")

            nalleles0 = np.count_nonzero(
                ~genotypedata_RR.loc[
                    genotypedata_RR["Sample.ID"].str.contains(f"{ID} Day 0"),
                    locicolumns,
                ].isna()
            )
            nallelesf = np.count_nonzero(
                ~genotypedata_RR.loc[
                    genotypedata_RR["Sample.ID"].str.contains(f"{ID} Day Failure"),
                    locicolumns,
                ].isna()
            )

            MOI0[i] = np.max([MOI0[i], nalleles0])
            MOIf[i] = np.max([MOIf[i], nallelesf])

    assert np.array_equal(MOI0, expected_MOI0), f"{MOI0} (expected {expected_MOI0})"
    assert np.array_equal(MOIf, expected_MOIf), f"{MOIf} (expected {expected_MOIf})"


def test_create_initial_state():
    # TODO: Finish implementing so test will pass
    expected_alleles0_firstCol = np.array([223.4, 229.3, 221.3, 261.7, 261.7, 238.4])
    expected_allelesf_firstCol = np.array([225.5, 243.6, 231.5, 239.6, 245.3, 219.4])
    expected_recoded0_firstCol = np.array([5, 12, 2, 4, 4, 8])
    expected_recodedf_firstCol = np.array([6, 3, 7, 9, 13, 1])

    maxMOI = 5
    ids = np.unique(
        ["BQ17-269_", "BD17-040_", "BD17-083_", "BD17-085_", "BD17-087_", "BD17-090_"]
    )
    locinames = np.unique(["X313", "X383", "TA1", "POLYA", "PFPK2", "X2490", "TA109"])
    alleles0 = np.zeros((ids.size, maxMOI * locinames.size))
    recoded0 = np.zeros((ids.size, maxMOI * locinames.size))
    allelesf = np.zeros((ids.size, maxMOI * locinames.size))
    recodedf = np.zeros((ids.size, maxMOI * locinames.size))

    for i, locus in enumerate(locinames):
        # locicolumns code is duplicated from MOI calculations
        locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

        oldalleles = genotypedata_RR.loc[:, locicolumns].to_numpy()
        """
        # TODO: What is this code doing?
        if (len(oldalleles.shape[1]) == 0) {
            oldalleles = matrix(oldalleles,length(oldalleles),1)
        }
        """
        newalleles = np.copy(oldalleles)
        ncolumns = oldalleles.shape[1]
        """
        for j in range(ncolumns):
            # TODO: Can't test until I have a recodeallele implementation
            newalleles[:,j] = np.array(list(map(
                range(0, oldalleles.shape[0]),
                lambda x: recodeallele(alleles_definitions_RR[i], oldalleles[x,j]))))
        """
        newalleles[np.isnan(newalleles)] = 0
        oldalleles[np.isnan(oldalleles)] = 0

        oldalleles[newalleles == 0] = 0

        startColumn = maxMOI * i
        # TODO: Changed start/end for indexing reasons in Python vs R; double-check that's valid
        endColumnOldAllele = maxMOI * i + oldalleles.shape[1]
        endColumnNewAllele = maxMOI * i + newalleles.shape[1]
        alleles0[:, startColumn:endColumnOldAllele] = oldalleles[
            genotypedata_RR["Sample.ID"].str.contains("Day 0"), :
        ]
        allelesf[:, startColumn:endColumnOldAllele] = oldalleles[
            genotypedata_RR["Sample.ID"].str.contains("Day Failure"), :
        ]
        recoded0[:, startColumn:endColumnNewAllele] = newalleles[
            genotypedata_RR["Sample.ID"].str.contains("Day 0"), :
        ]
        recodedf[:, startColumn:endColumnNewAllele] = newalleles[
            genotypedata_RR["Sample.ID"].str.contains("Day Failure"), :
        ]

    assert np.array_equal(
        alleles0[:, 0], expected_alleles0_firstCol
    ), f"{alleles0[:,0]} (expected {expected_alleles0_firstCol})"
    assert np.array_equal(
        recoded0[:, 0], expected_recoded0_firstCol
    ), f"{recoded0[:,0]} (expected {expected_recoded0_firstCol})"
    assert np.array_equal(
        allelesf[:, 0], expected_allelesf_firstCol
    ), f"{allelesf[:,0]} (expected {expected_allelesf_firstCol})"
    assert np.array_equal(
        recodedf[:, 0], expected_recodedf_firstCol
    ), f"{recodedf[:,0]} (expected {expected_recodedf_firstCol})"


def test_recode_additional_neutral():
    expected_first_column = np.array([6, 10, 1, 3, 7, 5, 9, 4, 3, 1, 0, 2, 11, 2])

    maxMOI = 5
    locinames = np.unique(["X313", "X383", "TA1", "POLYA", "PFPK2", "X2490", "TA109"])
    # TODO: additional_neutral currently just stubbed (not using actual values)
    additional_neutral = np.zeros((14, 25))

    recoded_additional_neutral = np.zeros((14, maxMOI * locinames.size))
    # TODO: This is almost the exact same code as in create_initial_state (should refactor both into a common function)

    for i, locus in enumerate(locinames):
        locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

        oldalleles = additional_neutral[:, locicolumns]  # TODO: stub
        """
        # TODO: What is this code doing?
        if (len(oldalleles.shape[1]) == 0) {
            oldalleles = matrix(oldalleles,length(oldalleles),1)
        }
        """
        newalleles = np.copy(oldalleles)
        ncolumns = oldalleles.shape[1]
        """
        for j in range(ncolumns):
            # TODO: Can't test until I have a recodeallele implementation
            newalleles[:,j] = np.array(list(map(
                range(0, oldalleles.shape[0]),
                lambda x: recodeallele(alleles_definitions_RR[i], oldalleles[x,j]))))
        """
        newalleles[np.isnan(newalleles)] = 0
        oldalleles[np.isnan(oldalleles)] = 0

        oldalleles[newalleles == 0] = 0

        startColumn = maxMOI * i
        # TODO: Changed start/end for indexing reasons in Python vs R; double-check that's valid
        endColumnOldAllele = maxMOI * i + oldalleles.shape[1]
        recoded_additional_neutral[:, startColumn:endColumnOldAllele] = newalleles

    assert np.array_equal(
        recoded_additional_neutral[:, 0], expected_first_column
    ), f"{recoded_additional_neutral[:,0]} (expected {expected_first_column})"


def test_initialize_hidden_alleles():
    # TODO: How to do expected values for stochastic functions?
    # ===============================
    # Initialize
    maxMOI = 5
    nids = 6
    nloci = 7

    MOI0 = np.array([2, 3, 2, 3, 2, 3])
    MOIf = np.array([2, 2, 2, 3, 2, 2])

    alleles0 = np.zeros((nids, maxMOI * nloci))
    recoded0 = np.zeros((nids, maxMOI * nloci))
    hidden0 = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    allelesf = np.zeros((nids, maxMOI * nloci))
    recodedf = np.zeros((nids, maxMOI * nloci))
    hiddenf = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)

    # TODO: Not sure what the proper shape of these are (frequencies_RR seems to have different-sized rows, with the 1st row all ints and the 2nd all floats???)
    # TODO: These are stubbed (waiting on correct implementation)
    frequencies_RR = np.random.random_sample((3, 7, 19))
    frequencies_RR[0, :, 0] *= frequencies_RR.shape[1]
    alleles_definitions_RR = np.zeros((13, 7, 2))
    # ===============================

    # TODO: Figure out what this code is overall trying to do (replace non-0 elements with random values? Why is it important that the values are assigned from frequencies_RR?)
    for i in range(nids):
        for j in range(nloci):
            # TODO: Code almost duplicated between top/bottom portions; refactor into single function? (needs 9 inputs: maxMOI, nids, nloci, MOIarray, alleles/recoded/hidden array, alleles_definitions_RR, frequencies_RR)
            # TODO: Start/end of what?
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

            if nalleles0 > 0:
                hidden0[i, whichnotmissing0] = 0
            if nmissing0 > 0:
                newhiddenalleles0 = np.random.choice(
                    np.arange(
                        0, int(frequencies_RR[0, j, 0])
                    ),  # Select from first row (count of how many probabilities they are)
                    size=nmissing0,
                    replace=True,
                    p=frequencies_RR[1, j, 0 : int(frequencies_RR[0, j, 0])]
                    / frequencies_RR[1, j, 0 : int(frequencies_RR[0, j, 0])].sum(),
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
                        0, int(frequencies_RR[0, j, 0])
                    ),  # Select from first row (count of how many probabilities they are)
                    size=nmissingf,
                    replace=True,
                    p=frequencies_RR[1, j, 0 : int(frequencies_RR[0, j, 0])]
                    / frequencies_RR[1, j, 0 : int(frequencies_RR[0, j, 0])].sum(),
                )  # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
                recodedf[i, whichmissingf] = newhiddenallelesf
                # calculate row means
                allelesf[i, whichmissingf] = np.mean(alleles_definitions_RR[j], axis=1)[
                    newhiddenallelesf
                ]  # hidden alleles get mean allele length
                hiddenf[i, whichmissingf] = 1

        # TODO: Figure out how to actually test this w/ random sampling?
        assert False, "test_initialize_hidden_alleles() using stubbed data for now"


def test_create_dvect():
    # TODO: This is a stub, not the actual data
    nloci = 7
    alleles_definitions_RR = np.array(
        [
            [217, 263],
            [85, 171],
            [71.5, 203.5],
            [102.5, 180.5],
            [136.5, 202.5],
            [71.5, 89.5],
            [146.5, 197.5],
        ]
    )

    # TODO: What does dvect stand for?
    dvect = np.zeros(
        1
        + int(
            np.rint(
                # Get the range (max-min) of the first "nloci" rows, then the max of all those
                np.ptp(alleles_definitions_RR[0:nloci], axis=1).max()
            )
        )
    )
    dvect[1] = 0.75
    dvect[2] = 0.2
    dvect[3] = 0.05

    assert dvect.size == 133, f"Dvect size {dvect.size} (expected {133})"


def test_initialize_recrudesences():
    maxMOI = 5
    nids = 6
    nloci = 7

    MOI0 = np.array([2, 3, 2, 3, 2, 3])
    MOIf = np.array([2, 2, 2, 3, 2, 2])

    recr0 = np.full_like(np.empty((nids, nloci)), np.nan)
    recrf = np.full_like(np.empty((nids, nloci)), np.nan)
    recr_repeats0 = np.full_like(np.empty((nids, nloci)), np.nan)
    recr_repeatsf = np.full_like(np.empty((nids, nloci)), np.nan)

    mindistance = np.zeros((nids, nloci))
    alldistance = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), np.nan)
    allrecrf = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), np.nan)
    classification = np.repeat(0, nids)

    # TODO: Stubbed data
    alleles0 = 100 * np.random.random_sample((nids, maxMOI * nloci))
    allelesf = 100 * np.random.random_sample((nids, maxMOI * nloci))
    recoded0 = np.random.randint(10, size=(nids, maxMOI * nloci))
    recodedf = np.random.randint(10, size=(nids, maxMOI * nloci))

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
                recoded0[i, maxMOI * j : maxMOI * (j + 1)]
                == recoded0[i, int(recr0[i, j])]
            )  # TODO: int() only needed for stub
            recr_repeatsf[i, j] = np.sum(
                recodedf[i, maxMOI * j : maxMOI * (j + 1)]
                == recodedf[i, int(recrf[i, j])]
            )  # TODO: int() only needed for stub
    # TODO: Actually test this


def test_correction_factor():
    expected_fist_row_col = np.array([0, 2, 24, 42, 4, 6, 12, 18, 20, 30, 2, 10, 26])

    nloci = 7
    # TODO: Confirm what the right alleles_definitions_RR shape is?
    alleles_definitions_RR = np.zeros((nloci, 13, 2))
    alleles_definitions_RR[0, :, :] = np.array(
        [
            [219, 221],
            [221, 223],
            [243, 245],
            [261, 263],
            [223, 225],
            [225, 227],
            [231, 233],
            [237, 239],
            [239, 241],
            [249, 251],
            [217, 219],
            [229, 231],
            [245, 247],
        ]
    )

    correction_distance_matrix = np.zeros(
        (nloci, alleles_definitions_RR.shape[1], alleles_definitions_RR.shape[1])
    )  # for each locus, matrix of distances between each allele
    for i in range(nloci):
        # Wrap mean call in "array" so we get a 2D array we can transpose (getting us a grid of distances, not just a 1D vector)
        distances = np.array([np.mean(alleles_definitions_RR[i], axis=1)])
        distance_combinations = np.abs(distances.T - distances)
        correction_distance_matrix[i] = distance_combinations

    assert np.array_equal(
        correction_distance_matrix.shape, np.array([7, 13, 13])
    ), f"{correction_distance_matrix.shape} (expected {np.array([7, 13, 13])})"
    assert np.array_equal(
        correction_distance_matrix[0, :, 0], expected_fist_row_col
    ), f"Row {correction_distance_matrix[0, :, 0]} (expected {expected_fist_row_col})"
    assert np.array_equal(
        correction_distance_matrix[0, 0, :], expected_fist_row_col
    ), f"Column {correction_distance_matrix[0, 0, :]} (expected {expected_fist_row_col})"


def test_run_mcmc_new_proposal():
    # TODO: Find what the actual expected ratio is for the stubbed inputs?
    expected_likelihood_ratio = np.zeros(6)

    maxMOI = 5
    nids = 6
    nloci = 7

    # TODO: Stubbed data
    frequencies_RR = np.random.random_sample((3, 7, 19))
    frequencies_RR[0, :, 0] = np.array([13, 16, 11, 13, 19, 5, 13])

    alldistance = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), 1)
    allrecrf = np.full_like(np.empty((nids, nloci, maxMOI ** 2)), 3).astype(int)
    hidden0 = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    hiddenf = np.full_like(np.empty((nids, maxMOI * nloci)), np.nan)
    classification = np.repeat(0, nids)

    dvect = np.ones(133)  # TODO: What to do when dvect has a 0 value?
    # TODO: What size is the correction matrix??? Appears to be different for each locus??? (set to max value in frequencies_rr)
    correction_distance_matrix = np.zeros((nloci, 19, 19)).astype(int)

    # propose new classification
    likelihoodratio = np.zeros(nids)
    # TODO: Finish vectorizing this
    for x in range(nids):
        # id mean for what?
        id_means = np.zeros(nloci)
        for y in range(nloci):
            id_means[y] = np.nanmean(
                dvect[np.round(alldistance[x, y, :]).astype(int)]
                # Should get an array of maxMOI**2 sums
                / np.sum(
                    frequencies_RR[1, y, : frequencies_RR.astype(int)[0, y, 0]]
                    * dvect[
                        correction_distance_matrix[  # TODO: Make sure multiplications are down the right axis (I believe they default to down columns, which is what I want)
                            y,
                            : frequencies_RR.astype(int)[
                                0, y, 0
                            ],  # TODO: Figure out how to do this more cleanly? (R impementation just used ":", assumed array had correct dimensions)
                            allrecrf[x, y, : maxMOI ** 2].astype(int),
                        ]
                    ],
                    axis=1,
                ),  # TODO: Verify it's the right axis?
            )
        likelihoodratio[x] = np.exp(np.sum(np.log(id_means)))

    z = np.random.uniform(size=nids)
    newclassification = classification
    newclassification[np.logical_and(classification == 0, z < likelihoodratio)] = 1
    newclassification[np.logical_and(classification == 1, z < 1 / likelihoodratio)] = 0
    classification = newclassification

    # propose new hidden states
    """
    # TODO: What does switch_hidden do? Is it entirely side effects? (Also,: can't run this yet, still waiting on implementation)
    for i in range(nids):
        switch_hidden(i)
    """

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

def test_update_dvect():
    nids = 6
    nloci = 7

    classification = np.repeat(0, nids)
    classification[0] = 1
    mindistance = np.zeros((nids, nloci))

    dvect = np.ones(133)

    # Inputs needed: classification, mindistance, dvect
    if np.sum(classification==1) >= 1:
        d_prior_alpha = 0
        d_prior_beta = 0
        d_posterior_alpha = d_prior_alpha + mindistance[classification==1,:].size
        d_posterior_beta = d_prior_beta + np.sum(np.round(mindistance[classification==1,:]))
        if (d_posterior_beta == 0):
            d_posterior_beta = np.sum(mindistance[classification==1,:])
        if (d_posterior_beta == 0): ## algorithm will get stuck if dposterior is allowed to go to 1 (TODO: Wait, so why is it setting d_posterior_beta to 1??)
            d_posterior_beta = 1

        dposterior = np.random.beta(d_posterior_alpha, d_posterior_beta)
        dvect = dposterior * (np.array(1-dposterior)**np.arange(0, dvect.size))
        dvect = dvect / np.sum(dvect)

    # TODO: Assert that dvect is updated correctly w/ beta function?

# =============================================================================

np.random.seed(0)

test_max_MOI()
test_getting_ids()
test_getting_locinames()
test_calculate_MOI()
# test_create_initial_state()       # TODO: Code not fully implemented yet
# test_recode_additional_neutral()  # TODO: Code not fully implemented yet
# test_initialize_hidden_alleles()  # TODO: Code not fully implemented yet
test_create_dvect()
test_initialize_recrudesences()
test_correction_factor()
test_run_mcmc_new_proposal()
test_update_dvect()
