#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest

import af2_analysis
from af2_analysis import docking
from .data_files import TEST_FILE_PATH


def test_cf_1_5_5_relax():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af2_analysis.Data(data_path)

    docking.pdockq2_lig(my_data)
    
    precision = 0.001
    expected_pdockq2 = [
        0.01199,
        0.01062,
        0.01429,
        0.01636,
        0.16844,
        0.01661,
        0.0115,
        0.02388,
        0.00895,
        0.02964,
        0.01174,
        0.01256,
        0.01341,
        0.05058,
        0.20396,
        0.01398,
        0.01116,
        0.01292,
        0.02662,
        0.19487,
        0.01218,
        0.01122,
        0.01373,
        0.05064,
        0.18241,
        0.01205,
        0.01135,
        0.01352,
        0.01025,
        0.02231,
        0.01376,
        0.0106,
        0.01412,
        0.04422,
        0.19023,
        0.01247,
        0.01078,
        0.01258,
        0.02291,
        0.20544,
    ]
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq2_lig"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    docking.pae_pep(my_data)

    expected_pae_pep_rec = [21.0262, 22.7601, 19.2283, 19.3637, 11.1904, 19.4724, 22.6415, 16.9034, 22.753, 15.5385, 21.0008, 22.4334, 19.7371, 15.1822, 9.361, 21.1361, 22.7517, 19.6134, 17.1336, 9.8214, 21.3164, 22.6982, 19.5113, 15.1621, 10.7245, 21.5201, 22.7322, 19.5693, 22.8056, 15.8663, 20.9644, 22.6516, 19.3624, 15.7163, 9.9666, 21.2954, 22.6243, 19.5583, 17.7726, 9.9227]

    assert np.all(
        [
            my_data.df.iloc[i]["PAE_pep_rec"]
            == pytest.approx(expected_pae_pep_rec[i], precision)
            for i in range(len(my_data.df))
        ]
    )
    expected_pae_rec_pep = [20.9843, 22.7467, 19.258, 19.289, 11.0753, 19.4125, 22.6326, 16.8903, 22.8446, 15.4945, 21.0322, 22.4293, 19.7231, 15.1634, 9.4072, 21.1865, 22.7387, 19.5859, 17.1389, 9.7928, 21.2765, 22.7175, 19.5033, 15.1969, 10.6566, 21.4941, 22.7618, 19.5856, 22.8432, 15.9111, 20.9533, 22.6817, 19.3914, 15.7483, 10.0023, 21.2943, 22.6347, 19.5152, 17.8586, 9.9088]

    assert np.all(
        [
            my_data.df.iloc[i]["PAE_rec_pep"]
            == pytest.approx(expected_pae_rec_pep[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    docking.LIS_pep(my_data)

    expected_LIS_pep_rec = [0.11039, 0.10359, 0.08238, 0.10013, 0.41793, 0.16121, 0.11254, 0.17166, 0.02854, 0.19839, 0.10146, 0.12784, 0.08725, 0.19706, 0.46334, 0.12865, 0.09632, 0.09369, 0.13473, 0.44639, 0.11678, 0.11155, 0.09551, 0.18691, 0.43507, 0.11824, 0.09029, 0.09923, 0.02556, 0.18627, 0.12552, 0.10105, 0.09658, 0.17747, 0.44164, 0.11034, 0.10487, 0.10166, 0.12171, 0.48107]
    # print([round(val, 5) for val in my_data.df["LIS_pep_rec"]])

    assert np.all(
        [
            my_data.df.iloc[i]["LIS_pep_rec"]
            == pytest.approx(expected_LIS_pep_rec[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    print([round(val, 5) for val in my_data.df["LIS_rec_pep"]])
    expected_LIS_rec_pep = [0.10748, 0.10333, 0.08494, 0.0996, 0.41781, 0.15442, 0.11433, 0.17288, 0.03479, 0.20141, 0.10314, 0.12448, 0.08833, 0.1976, 0.45977, 0.1319, 0.09803, 0.09379, 0.13357, 0.44645, 0.11865, 0.11202, 0.09505, 0.18964, 0.43426, 0.11507, 0.09413, 0.10194, 0.03833, 0.18998, 0.12596, 0.10077, 0.09478, 0.17737, 0.44065, 0.11144, 0.10867, 0.10031, 0.12168, 0.48279]
    assert np.all(
        [
            my_data.df.iloc[i]["LIS_rec_pep"]
            == pytest.approx(expected_LIS_rec_pep[i], precision)
            for i in range(len(my_data.df))
        ]
    )