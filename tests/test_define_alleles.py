import os
import sys

import pandas as pd
import pytest
import numpy as np

import api.define_alleles as Define_Alleles
from api.algorithm_instance import AlgorithmInstance
import api.recrudescence_utils as recrudescence_utils
from api.recrudescence_file_parser import RecrudescenceFileParser

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

locirepeats = np.array([2, 2, 3, 3, 3, 3, 3])
maxk = np.array([30, 30, 30, 30, 30, 30, 30])

@pytest.fixture
def expected_Benguela_define_alleles():
    expected_list = []
    expected_columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[219.0, 221.0],[221.0, 223.0],[243.0, 245.0],[261.0, 263.0],[223.0, 225.0],
         [225.0, 227.0],[231.0, 233.0],[237.0, 239.0],[239.0, 241.0],[249.0, 251.0],
         [217.0, 219.0],[229.0, 231.0],[245.0, 247.0]],
         columns=expected_columns)
    index_1 = pd.DataFrame(
        [[123.0, 125.0],[139.0, 141.0],[103.0, 105.0],[85.0, 87.0],[121.0, 123.0],
         [143.0, 145.0],[145.0, 147.0],[87.0, 89.0],[135.0, 137.0],[137.0, 139.0],
         [147.0, 149.0],[149.0, 151.0],[151.0, 153.0],[161.0, 163.0],[163.0, 165.0],
         [169.0, 171.0]],
         columns=expected_columns)
    index_2 = pd.DataFrame(
        [[176.5, 179.5],[161.5, 164.5],[164.5, 167.5],[167.5, 170.5],[173.5, 176.5],
         [170.5, 173.5],[71.5, 74.5],[179.5, 182.5],[158.5, 161.5],[182.5, 185.5],
         [200.5, 203.5]],
        columns=expected_columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[150.5, 153.5],[162.5, 165.5],[165.5, 168.5],[156.5, 159.5],
         [141.5, 144.5],[168.5, 171.5],[171.5, 174.5],[102.5, 105.5],[111.5, 114.5],
         [159.5, 162.5],[174.5, 177.5],[177.5, 180.5]],
         columns=expected_columns)
    index_4 = pd.DataFrame(
        [[160.5, 163.5],[166.5, 169.5],[169.5, 172.5],[175.5, 178.5],[172.5, 175.5],
         [163.5, 166.5],[178.5, 181.5],[148.5, 151.5],[154.5, 157.5],[157.5, 160.5],
         [184.5, 187.5],[193.5, 196.5],[136.5, 139.5],[139.5, 142.5],[181.5, 184.5],
         [187.5, 190.5],[190.5, 193.5],[196.5, 199.5],[199.5, 202.5]],
         columns=expected_columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[77.5, 80.5],[71.5, 74.5],[83.5, 86.5],[86.5, 89.5]],
        columns=expected_columns)
    index_6 = pd.DataFrame(
        [[161.5, 164.5],[158.5, 161.5],[146.5, 149.5],[173.5, 176.5],[170.5, 173.5],
         [149.5, 152.5],[167.5, 170.5],[164.5, 167.5],[176.5, 179.5],[179.5, 182.5],
         [182.5, 185.5],[185.5, 188.5],[194.5, 197.5]],
         columns=expected_columns)

    expected_list.append(index_0)
    expected_list.append(index_1)
    expected_list.append(index_2)
    expected_list.append(index_3)
    expected_list.append(index_4)
    expected_list.append(index_5)
    expected_list.append(index_6)

    return expected_list

@pytest.fixture
def expected_Lunda_Sul_define_alleles():
    expected_list = []
    expected_columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[218.0, 220.0],[244.0, 246.0],[258.0, 260.0],[216.0, 218.0],[224.0, 226.0],
         [232.0, 234.0],[234.0, 236.0],[236.0, 238.0],[210.0, 212.0],[220.0, 222.0],
         [230.0, 232.0],[238.0, 240.0],[240.0, 242.0],[246.0, 248.0],[248.0, 250.0],
         [256.0, 258.0],[274.0, 276.0]],
         columns=expected_columns)
    index_1 = pd.DataFrame(
        [[123.0, 125.0],[139.0, 141.0],[85.0, 87.0],[137.0, 139.0],[143.0, 145.0],
         [101.0, 103.0],[121.0, 123.0],[135.0, 137.0],[145.0, 147.0],[87.0, 89.0],
         [99.0, 101.0],[107.0, 109.0],[127.0, 129.0],[149.0, 151.0],[151.0, 153.0],
         [153.0, 155.0],[155.0, 157.0],[169.0, 171.0]],
         columns=expected_columns)
    index_2 = pd.DataFrame(
        [[158.5, 161.5],[164.5, 167.5],[167.5, 170.5],[170.5, 173.5],[173.5, 176.5],
         [176.5, 179.5],[179.5, 182.5],[71.5, 74.5],[161.5, 164.5],[185.5, 188.5]],
        columns=expected_columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[162.5, 165.5],[156.5, 159.5],[147.5, 150.5],[150.5, 153.5],
         [168.5, 171.5],[174.5, 177.5],[159.5, 162.5],[165.5, 168.5],[135.5, 138.5],
         [141.5, 144.5],[180.5, 183.5],[186.5, 189.5]],
         columns=expected_columns)
    index_4 = pd.DataFrame(
        [[163.5, 166.5],[160.5, 163.5],[166.5, 169.5],[157.5, 160.5],[169.5, 172.5],
         [175.5, 178.5],[172.5, 175.5],[154.5, 157.5],[178.5, 181.5],[184.5, 187.5],
         [181.5, 184.5],[187.5, 190.5],[190.5, 193.5]],
         columns=expected_columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[71.5, 74.5],[77.5, 80.5]],
        columns=expected_columns)
    index_6 = pd.DataFrame(
        [[161.5, 164.5],[158.5, 161.5],[173.5, 176.5],[170.5, 173.5],[146.5, 149.5],
         [149.5, 152.5],[164.5, 167.5],[152.5, 155.5],[179.5, 182.5],[167.5, 170.5],
         [176.5, 179.5],[200.5, 203.5],[203.5, 206.5],[206.5, 209.5],[209.5, 212.5]],
         columns=expected_columns)

    expected_list.append(index_0)
    expected_list.append(index_1)
    expected_list.append(index_2)
    expected_list.append(index_3)
    expected_list.append(index_4)
    expected_list.append(index_5)
    expected_list.append(index_6)

    return expected_list

@pytest.fixture
def expected_Zaire_define_alleles():
    expected_list = []
    expected_columns = ['0', '1']
    index_0 = pd.DataFrame(
        [[216.0, 218.0],[224.0, 226.0],[232.0, 234.0],[220.0, 222.0],[228.0, 230.0],
         [242.0, 244.0],[238.0, 240.0],[240.0, 242.0],[222.0, 224.0],[230.0, 232.0],
         [250.0, 252.0],[210.0, 212.0],[214.0, 216.0],[234.0, 236.0],[236.0, 238.0],
         [244.0, 246.0],[254.0, 256.0],[260.0, 262.0]],
         columns=expected_columns)
    index_1 = pd.DataFrame(
        [[137.0, 139.0],[139.0, 141.0],[123.0, 125.0],[101.0, 103.0],[135.0, 137.0],
         [83.0, 85.0],[121.0, 123.0],[141.0, 143.0],[85.0, 87.0],[97.0, 99.0],
         [99.0, 101.0],[129.0, 131.0],[163.0, 165.0],[171.0, 173.0],[87.0, 89.0],
         [89.0, 91.0],[103.0, 105.0],[105.0, 107.0],[143.0, 145.0],[147.0, 149.0],
         [149.0, 151.0],[151.0, 153.0],[161.0, 163.0]],
         columns=expected_columns)
    index_2 = pd.DataFrame(
        [[158.5, 161.5],[164.5, 167.5],[170.5, 173.5],[176.5, 179.5],[161.5, 164.5],
         [167.5, 170.5],[182.5, 185.5],[137.5, 140.5],[173.5, 176.5],[179.5, 182.5],
         [134.5, 137.5],[140.5, 143.5],[185.5, 188.5],[191.5, 194.5]],
        columns=expected_columns)
    index_3 = pd.DataFrame(
        [[153.5, 156.5],[150.5, 153.5],[165.5, 168.5],[147.5, 150.5],[177.5, 180.5],
         [171.5, 174.5],[141.5, 144.5],[156.5, 159.5],[159.5, 162.5],[138.5, 141.5],
         [162.5, 165.5],[168.5, 171.5],[117.5, 120.5],[135.5, 138.5],[180.5, 183.5]],
         columns=expected_columns)
    index_4 = pd.DataFrame(
        [[160.5, 163.5],[166.5, 169.5],[169.5, 172.5],[157.5, 160.5],[163.5, 166.5],
         [172.5, 175.5],[175.5, 178.5],[178.5, 181.5],[181.5, 184.5],[151.5, 154.5],
         [184.5, 187.5],[187.5, 190.5],[133.5, 136.5],[136.5, 139.5],[142.5, 145.5],
         [145.5, 148.5]],
         columns=expected_columns)
    index_5 = pd.DataFrame(
        [[80.5, 83.5],[77.5, 80.5],[83.5, 86.5],[71.5, 74.5],],
        columns=expected_columns)
    index_6 = pd.DataFrame(
        [[170.5, 173.5],[158.5, 161.5],[161.5, 164.5],[173.5, 176.5],[149.5, 152.5],
         [182.5, 185.5],[146.5, 149.5],[164.5, 167.5],[176.5, 179.5],[185.5, 188.5],
         [143.5, 146.5],[152.5, 155.5],[155.5, 158.5],[167.5, 170.5],[179.5, 182.5],
         [194.5, 197.5]],
         columns=expected_columns)

    expected_list.append(index_0)
    expected_list.append(index_1)
    expected_list.append(index_2)
    expected_list.append(index_3)
    expected_list.append(index_4)
    expected_list.append(index_5)
    expected_list.append(index_6)

    return expected_list

def test_Benguela__define_alleles_output_correct(expected_Benguela_define_alleles):
    result_list = Define_Alleles.define_alleles(
        pd.concat([genotypedata_RR_Benguela, additional_neutral_Benguela]), locirepeats, maxk
    )

    # assert result_list.length == expected_Benguela_define_alleles.length

    assert len(result_list) == len(expected_Benguela_define_alleles)
    
    for index in range(len(result_list)):
        pd.testing.assert_frame_equal(result_list[index], expected_Benguela_define_alleles[index])

def test_Lunda_Sul__define_alleles_output_correct(expected_Lunda_Sul_define_alleles):
    result_list = Define_Alleles.define_alleles(
        pd.concat([genotypedata_RR_Lunda_Sul, additional_neutral_Lunda_Sul]), locirepeats, maxk
    )

    # assert result_list.length == expected_Benguela_define_alleles.length

    assert len(result_list) == len(expected_Lunda_Sul_define_alleles)
    
    for index in range(len(result_list)):
        pd.testing.assert_frame_equal(result_list[index], expected_Lunda_Sul_define_alleles[index])

def test_Zaire__define_alleles_output_correct(expected_Zaire_define_alleles):
    result_list = Define_Alleles.define_alleles(
        pd.concat([genotypedata_RR_Zaire, additional_neutral_Zaire]), locirepeats, maxk
    )

    # assert result_list.length == expected_Benguela_define_alleles.length

    assert len(result_list) == len(expected_Zaire_define_alleles)
    
    for index in range(len(result_list)):
        pd.testing.assert_frame_equal(result_list[index], expected_Zaire_define_alleles[index])