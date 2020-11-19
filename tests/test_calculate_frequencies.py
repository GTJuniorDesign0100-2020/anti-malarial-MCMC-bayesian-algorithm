import os
import sys

import pandas as pd
import pytest
import numpy as np

from api.algorithm_instance import AlgorithmInstance
import api.recrudescence_utils as recrudescence_utils
from api.recrudescence_file_parser import RecrudescenceFileParser
from api.calculate_frequencies import calculate_frequencies3, Frequencies

# NOTE: Makes this reliant on AlgorithmInstance tests passing
example_file = os.path.join(
    os.path.dirname(__file__),
    '../Angola2017_example.xlsx')
genotypedata, additional = RecrudescenceFileParser.parse_file(example_file)

genotypedata_RR_Benguela = AlgorithmInstance._get_samples_from_site(genotypedata, 'Benguela')
additional_neutral_Benguela = AlgorithmInstance._replace_sample_names(
    AlgorithmInstance._get_samples_from_site(additional, 'Benguela'),'Additional_')

genotypedata_RR_Lunda_Sul = AlgorithmInstance._get_samples_from_site(genotypedata, 'Lunda Sul')
additional_neutral_Lunda_Sul = AlgorithmInstance._replace_sample_names(
    AlgorithmInstance._get_samples_from_site(additional, 'Lunda Sul'),'Additional_')

genotypedata_RR_Zaire = AlgorithmInstance._get_samples_from_site(genotypedata, 'Zaire')
additional_neutral_Zaire = AlgorithmInstance._replace_sample_names(
    AlgorithmInstance._get_samples_from_site(additional, 'Zaire'),'Additional_')

@pytest.fixture
def Benguela_define_alleles():
    define_alleles = []
    columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[219.0, 221.0],[221.0, 223.0],[243.0, 245.0],[261.0, 263.0],[223.0, 225.0],
         [225.0, 227.0],[231.0, 233.0],[237.0, 239.0],[239.0, 241.0],[249.0, 251.0],
         [217.0, 219.0],[229.0, 231.0],[245.0, 247.0]],
         columns=columns)
    index_1 = pd.DataFrame(
        [[123.0, 125.0],[139.0, 141.0],[103.0, 105.0],[85.0, 87.0],[121.0, 123.0],
         [143.0, 145.0],[145.0, 147.0],[87.0, 89.0],[135.0, 137.0],[137.0, 139.0],
         [147.0, 149.0],[149.0, 151.0],[151.0, 153.0],[161.0, 163.0],[163.0, 165.0],
         [169.0, 171.0]],
         columns=columns)
    index_2 = pd.DataFrame(
        [[176.5, 179.5],[161.5, 164.5],[164.5, 167.5],[167.5, 170.5],[173.5, 176.5],
         [170.5, 173.5],[71.5, 74.5],[179.5, 182.5],[158.5, 161.5],[182.5, 185.5],
         [200.5, 203.5]],
        columns=columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[150.5, 153.5],[162.5, 165.5],[165.5, 168.5],[156.5, 159.5],
         [141.5, 144.5],[168.5, 171.5],[171.5, 174.5],[102.5, 105.5],[111.5, 114.5],
         [159.5, 162.5],[174.5, 177.5],[177.5, 180.5]],
         columns=columns)
    index_4 = pd.DataFrame(
        [[160.5, 163.5],[166.5, 169.5],[169.5, 172.5],[175.5, 178.5],[172.5, 175.5],
         [163.5, 166.5],[178.5, 181.5],[148.5, 151.5],[154.5, 157.5],[157.5, 160.5],
         [184.5, 187.5],[193.5, 196.5],[136.5, 139.5],[139.5, 142.5],[181.5, 184.5],
         [187.5, 190.5],[190.5, 193.5],[196.5, 199.5],[199.5, 202.5]],
         columns=columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[77.5, 80.5],[71.5, 74.5],[83.5, 86.5],[86.5, 89.5]],
        columns=columns)
    index_6 = pd.DataFrame(
        [[161.5, 164.5],[158.5, 161.5],[146.5, 149.5],[173.5, 176.5],[170.5, 173.5],
         [149.5, 152.5],[167.5, 170.5],[164.5, 167.5],[176.5, 179.5],[179.5, 182.5],
         [182.5, 185.5],[185.5, 188.5],[194.5, 197.5]],
         columns=columns)

    define_alleles.append(index_0)
    define_alleles.append(index_1)
    define_alleles.append(index_2)
    define_alleles.append(index_3)
    define_alleles.append(index_4)
    define_alleles.append(index_5)
    define_alleles.append(index_6)
    return define_alleles

@pytest.fixture
def Lunda_Sul_define_alleles():
    define_alleles = []
    columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[218.0, 220.0],[244.0, 246.0],[258.0, 260.0],[216.0, 218.0],[224.0, 226.0],
         [232.0, 234.0],[234.0, 236.0],[236.0, 238.0],[210.0, 212.0],[220.0, 222.0],
         [230.0, 232.0],[238.0, 240.0],[240.0, 242.0],[246.0, 248.0],[248.0, 250.0],
         [256.0, 258.0],[274.0, 276.0]],
         columns=columns)
    index_1 = pd.DataFrame(
        [[123.0, 125.0],[139.0, 141.0],[85.0, 87.0],[137.0, 139.0],[143.0, 145.0],
         [101.0, 103.0],[121.0, 123.0],[135.0, 137.0],[145.0, 147.0],[87.0, 89.0],
         [99.0, 101.0],[107.0, 109.0],[127.0, 129.0],[149.0, 151.0],[151.0, 153.0],
         [153.0, 155.0],[155.0, 157.0],[169.0, 171.0]],
         columns=columns)
    index_2 = pd.DataFrame(
        [[158.5, 161.5],[164.5, 167.5],[167.5, 170.5],[170.5, 173.5],[173.5, 176.5],
         [176.5, 179.5],[179.5, 182.5],[71.5, 74.5],[161.5, 164.5],[185.5, 188.5]],
        columns=columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[162.5, 165.5],[156.5, 159.5],[147.5, 150.5],[150.5, 153.5],
         [168.5, 171.5],[174.5, 177.5],[159.5, 162.5],[165.5, 168.5],[135.5, 138.5],
         [141.5, 144.5],[180.5, 183.5],[186.5, 189.5]],
         columns=columns)
    index_4 = pd.DataFrame(
        [[163.5, 166.5],[160.5, 163.5],[166.5, 169.5],[157.5, 160.5],[169.5, 172.5],
         [175.5, 178.5],[172.5, 175.5],[154.5, 157.5],[178.5, 181.5],[184.5, 187.5],
         [181.5, 184.5],[187.5, 190.5],[190.5, 193.5]],
         columns=columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[71.5, 74.5],[77.5, 80.5]],
        columns=columns)
    index_6 = pd.DataFrame(
        [[161.5, 164.5],[158.5, 161.5],[173.5, 176.5],[170.5, 173.5],[146.5, 149.5],
         [149.5, 152.5],[164.5, 167.5],[152.5, 155.5],[179.5, 182.5],[167.5, 170.5],
         [176.5, 179.5],[200.5, 203.5],[203.5, 206.5],[206.5, 209.5],[209.5, 212.5]],
         columns=columns)

    define_alleles.append(index_0)
    define_alleles.append(index_1)
    define_alleles.append(index_2)
    define_alleles.append(index_3)
    define_alleles.append(index_4)
    define_alleles.append(index_5)
    define_alleles.append(index_6)
    return define_alleles

@pytest.fixture
def Zaire_define_alleles():
    define_alleles = []
    columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[216.0, 218.0],[224.0, 226.0],[232.0, 234.0],[220.0, 222.0],[228.0, 230.0],
         [242.0, 244.0],[238.0, 240.0],[240.0, 242.0],[222.0, 224.0],[230.0, 232.0],
         [250.0, 252.0],[210.0, 212.0],[214.0, 216.0],[234.0, 236.0],[236.0, 238.0],
         [244.0, 246.0],[254.0, 256.0],[260.0, 262.0]],
         columns=columns)
    index_1 = pd.DataFrame(
        [[137.0, 139.0],[139.0, 141.0],[123.0, 125.0],[101.0, 103.0],[135.0, 137.0],
         [83.0, 85.0],[121.0, 123.0],[141.0, 143.0],[85.0, 87.0],[97.0, 99.0],
         [99.0, 101.0],[129.0, 131.0],[163.0, 165.0],[171.0, 173.0],[87.0, 89.0],
         [89.0, 91.0],[103.0, 105.0],[105.0, 107.0],[143.0, 145.0],[147.0, 149.0],
         [149.0, 151.0],[151.0, 153.0],[161.0, 163.0]],
         columns=columns)
    index_2 = pd.DataFrame(
        [[158.5, 161.5],[164.5, 167.5],[170.5, 173.5],[176.5, 179.5],[161.5, 164.5],
         [167.5, 170.5],[182.5, 185.5],[137.5, 140.5],[173.5, 176.5],[179.5, 182.5],
         [134.5, 137.5],[140.5, 143.5],[185.5, 188.5],[191.5, 194.5]],
        columns=columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[150.5, 153.5],[165.5, 168.5],[147.5, 150.5],[177.5, 180.5],
         [171.5, 174.5],[141.5, 144.5],[156.5, 159.5],[159.5, 162.5],[138.5, 141.5],
         [162.5, 165.5],[168.5, 171.5],[117.5, 120.5],[135.5, 138.5],[180.5, 183.5]],
         columns=columns)
    index_4 = pd.DataFrame(
        [[160.5, 163.5],[166.5, 169.5],[169.5, 172.5],[157.5, 160.5],[163.5, 166.5],
         [172.5, 175.5],[175.5, 178.5],[178.5, 181.5],[181.5, 184.5],[151.5, 154.5],
         [184.5, 187.5],[187.5, 190.5],[133.5, 136.5],[136.5, 139.5],[142.5, 145.5],
         [145.5, 148.5]],
         columns=columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[77.5, 80.5],[83.5, 86.5],[71.5, 74.5],],
        columns=columns)
    index_6 = pd.DataFrame(
        [[170.5, 173.5],[158.5, 161.5],[161.5, 164.5],[173.5, 176.5],[149.5, 152.5],
         [182.5, 185.5],[146.5, 149.5],[164.5, 167.5],[176.5, 179.5],[185.5, 188.5],
         [143.5, 146.5],[152.5, 155.5],[155.5, 158.5],[167.5, 170.5],[179.5, 182.5],
         [194.5, 197.5]],
         columns=columns)

    define_alleles.append(index_0)
    define_alleles.append(index_1)
    define_alleles.append(index_2)
    define_alleles.append(index_3)
    define_alleles.append(index_4)
    define_alleles.append(index_5)
    define_alleles.append(index_6)
    return define_alleles


@pytest.fixture
def expected_Benguela_frequencies():
    expected_list = []

    index_0 = np.array([13, 16, 11, 13, 19, 5, 13])
    index_1 = np.array([
        [0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.07407407,
        0.07407407, 0.07407407, 0.07407407, 0.07407407, 0.07407407,
        0.03703704, 0.03703704, 0.03703704, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
        [0.24324324, 0.21621622, 0.08108108, 0.05405405, 0.05405405,
        0.05405405, 0.05405405, 0.02702703, 0.02702703, 0.02702703,
        0.02702703, 0.02702703, 0.02702703, 0.02702703, 0.02702703,
        0.02702703, 0.        , 0.        , 0.        ],
        [0.2       , 0.17142857, 0.11428571, 0.11428571, 0.11428571,
        0.08571429, 0.05714286, 0.05714286, 0.02857143, 0.02857143,
        0.02857143, 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
        [0.18181818, 0.15151515, 0.12121212, 0.12121212, 0.09090909,
        0.06060606, 0.06060606, 0.06060606, 0.03030303, 0.03030303,
        0.03030303, 0.03030303, 0.03030303, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
        [0.1372549 , 0.11764706, 0.11764706, 0.09803922, 0.07843137,
        0.05882353, 0.05882353, 0.03921569, 0.03921569, 0.03921569,
        0.03921569, 0.03921569, 0.01960784, 0.01960784, 0.01960784,
        0.01960784, 0.01960784, 0.01960784, 0.01960784],
        [0.57142857, 0.25      , 0.10714286, 0.03571429, 0.03571429,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
        [0.23076923, 0.19230769, 0.13461538, 0.13461538, 0.09615385,
        0.05769231, 0.03846154, 0.01923077, 0.01923077, 0.01923077,
        0.01923077, 0.01923077, 0.01923077, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ]
        ])
    index_2 = np.array([0.4821891, 0.14313027, 0.08446277, 0.24941966, 0.30618052, 0.25014729, 0.79079775])

    expected_frequencies_RR = Frequencies(index_0, index_1, index_2)
    return expected_frequencies_RR

@pytest.fixture
def expected_Lunda_Sul_frequencies():
    expected_list = []

    index_0 = np.array([17, 18, 10, 13, 13,  3, 15])
    index_1 = np.array([
        [0.11764706, 0.11764706, 0.11764706, 0.08823529, 0.08823529,
        0.08823529, 0.05882353, 0.05882353, 0.02941176, 0.02941176,
        0.02941176, 0.02941176, 0.02941176, 0.02941176, 0.02941176,
        0.02941176, 0.02941176, 0.        ],
        [0.225     , 0.125     , 0.075     , 0.075     , 0.075     ,
        0.05      , 0.05      , 0.05      , 0.05      , 0.025     ,
        0.025     , 0.025     , 0.025     , 0.025     , 0.025     ,
        0.025     , 0.025     , 0.025     ],
        [0.22857143, 0.17142857, 0.11428571, 0.11428571, 0.11428571,
        0.11428571, 0.05714286, 0.02857143, 0.02857143, 0.02857143,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.20512821, 0.17948718, 0.1025641 , 0.07692308, 0.07692308,
        0.07692308, 0.07692308, 0.05128205, 0.05128205, 0.02564103,
        0.02564103, 0.02564103, 0.02564103, 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.2       , 0.13333333, 0.13333333, 0.08888889, 0.08888889,
        0.08888889, 0.06666667, 0.04444444, 0.04444444, 0.04444444,
        0.02222222, 0.02222222, 0.02222222, 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.66666667, 0.16666667, 0.16666667, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.22033898, 0.15254237, 0.15254237, 0.10169492, 0.06779661,
        0.06779661, 0.06779661, 0.03389831, 0.03389831, 0.01694915,
        0.01694915, 0.01694915, 0.01694915, 0.01694915, 0.01694915,
        0.        , 0.        , 0.        ]
        ])
    index_2 = np.array([0.30069772, 0.32102995, 0.35015972, 0.13598753, 0.18125558, 0.23458575, 0.50453622])

    expected_frequencies_RR = Frequencies(index_0, index_1, index_2)
    return expected_frequencies_RR

@pytest.fixture
def expected_Zaire_frequencies():
    expected_list = []

    index_0 = np.array([18, 23, 14, 15, 16,  4, 16])
    index_1 = np.array([
        [0.14545455, 0.12727273, 0.12727273, 0.10909091, 0.07272727,
        0.07272727, 0.05454545, 0.05454545, 0.03636364, 0.03636364,
        0.03636364, 0.01818182, 0.01818182, 0.01818182, 0.01818182,
        0.01818182, 0.01818182, 0.01818182, 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.14285714, 0.14285714, 0.12987013, 0.09090909, 0.09090909,
        0.05194805, 0.03896104, 0.03896104, 0.02597403, 0.02597403,
        0.02597403, 0.02597403, 0.02597403, 0.02597403, 0.01298701,
        0.01298701, 0.01298701, 0.01298701, 0.01298701, 0.01298701,
        0.01298701, 0.01298701, 0.01298701],
        [0.18032787, 0.13114754, 0.1147541 , 0.1147541 , 0.09836066,
        0.09836066, 0.06557377, 0.04918033, 0.04918033, 0.03278689,
        0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.19672131, 0.16393443, 0.1147541 , 0.08196721, 0.08196721,
        0.06557377, 0.04918033, 0.04918033, 0.04918033, 0.03278689,
        0.03278689, 0.03278689, 0.01639344, 0.01639344, 0.01639344,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.2375    , 0.1625    , 0.125     , 0.1       , 0.0875    ,
        0.05      , 0.0375    , 0.0375    , 0.0375    , 0.025     ,
        0.025     , 0.025     , 0.0125    , 0.0125    , 0.0125    ,
        0.0125    , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.62745098, 0.31372549, 0.03921569, 0.01960784, 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ],
        [0.19      , 0.18      , 0.15      , 0.12      , 0.07      ,
        0.07      , 0.06      , 0.04      , 0.04      , 0.02      ,
        0.01      , 0.01      , 0.01      , 0.01      , 0.01      ,
        0.01      , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        ]
        ])
    index_2 = np.array([0.22105065, 0.22734084, 0.24295935, 0.2933828 , 0.35520227,0.11494194, 0.47930888])

    expected_frequencies_RR = Frequencies(index_0, index_1, index_2)
    return expected_frequencies_RR


def test_Benguela_calculate_frequencies_output_correct(expected_Benguela_frequencies, Benguela_define_alleles):
    result_list = calculate_frequencies3(
        pd.concat([genotypedata_RR_Benguela, additional_neutral_Benguela]), Benguela_define_alleles)

    np.testing.assert_array_almost_equal(expected_Benguela_frequencies.lengths, result_list.lengths, decimal=5)
    np.testing.assert_array_almost_equal(expected_Benguela_frequencies.matrix, result_list.matrix, decimal=5)
    np.testing.assert_array_almost_equal(expected_Benguela_frequencies.variability, result_list.variability, decimal=5)

def test_Lunda_Sul_calculate_frequencies_output_correct(expected_Lunda_Sul_frequencies, Lunda_Sul_define_alleles):
    result_list = calculate_frequencies3(
        pd.concat([genotypedata_RR_Lunda_Sul, additional_neutral_Lunda_Sul]), Lunda_Sul_define_alleles)

    np.testing.assert_array_almost_equal(expected_Lunda_Sul_frequencies.lengths, result_list.lengths, decimal=5)
    np.testing.assert_array_almost_equal(expected_Lunda_Sul_frequencies.matrix, result_list.matrix, decimal=5)
    np.testing.assert_array_almost_equal(expected_Lunda_Sul_frequencies.variability, result_list.variability, decimal=5)

def test_Zaire_calculate_frequencies_output_correct(expected_Zaire_frequencies, Zaire_define_alleles):
    result_list = calculate_frequencies3(
        pd.concat([genotypedata_RR_Zaire, additional_neutral_Zaire]), Zaire_define_alleles)

    np.testing.assert_array_almost_equal(expected_Zaire_frequencies.lengths, result_list.lengths, decimal=5)
    np.testing.assert_array_almost_equal(expected_Zaire_frequencies.matrix, result_list.matrix, decimal=5)
    np.testing.assert_array_almost_equal(expected_Zaire_frequencies.variability, result_list.variability, decimal=5)

