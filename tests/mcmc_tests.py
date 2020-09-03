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

        startColumn = maxMOI * (
            i - 1
        )  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
        endColumnOldAllele = maxMOI * (i - 1) + oldalleles.shape[1]
        endColumnNewAllele = maxMOI * (i - 1) + newalleles.shape[1]
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

        startColumn = maxMOI * (
            i - 1
        )  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
        endColumnOldAllele = maxMOI * (i - 1) + oldalleles.shape[1]
        recoded_additional_neutral[:, startColumn:endColumnOldAllele] = newalleles

    assert np.array_equal(
        recoded_additional_neutral[:, 0], expected_first_column
    ), f"{recoded_additional_neutral[:,0]} (expected {expected_first_column})"


test_max_MOI()
test_getting_ids()
test_getting_locinames()
test_calculate_MOI()
# test_create_initial_state() # TODO: Code not fully implemented yet
# test_recode_additional_neutral() # TODO: Code not fully implemented yet
