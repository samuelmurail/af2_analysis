#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest

import af2_analysis
from af2_analysis import analysis
from .data_files import TEST_FILE_PATH


def test_cf_1_5_5_relax():
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af2_analysis.Data(data_path)
    analysis.pdockq(my_data)

    expected_pdockq = [
        0.03322,
        0.02051,
        0.03,
        0.02971,
        0.12094,
        0.02849,
        0.02246,
        0.04845,
        0.0198,
        0.07149,
        0.03319,
        0.02379,
        0.02763,
        0.05582,
        0.13833,
        0.02418,
        0.02112,
        0.02521,
        0.04189,
        0.14153,
        0.02952,
        0.02123,
        0.02851,
        0.05237,
        0.13703,
        0.02909,
        0.02042,
        0.02843,
        0.02069,
        0.08227,
        0.02312,
        0.02025,
        0.02824,
        0.05085,
        0.13916,
        0.0294,
        0.02063,
        0.02542,
        0.03623,
        0.14257,
    ]
    assert 2.2 == pytest.approx(2.3, 0.1)

    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq"] == pytest.approx(expected_pdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.pdockq2(my_data)

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
            my_data.df.iloc[i]["pdockq2_B"]
            == pytest.approx(expected_pdockq2[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    analysis.LIS_matrix(my_data)
    expected_LIS_0 = [[0.44785182, 0.11038826], [0.10748148, 0.44471074]]
    np.testing.assert_allclose(np.array(my_data.df["LIS"][0]), np.array(expected_LIS_0))

    analysis.inter_chain_pae(my_data)

    expected_PAE_A_B = [
        20.984252,
        22.746667,
        19.257993,
        19.288968,
        11.075300,
        19.412483,
        22.632647,
        16.890261,
        22.844563,
        15.494473,
        21.032171,
        22.429257,
        19.723124,
        15.163367,
        9.407222,
        21.186502,
        22.738668,
        19.585907,
        17.138900,
        9.792789,
        21.276542,
        22.717523,
        19.503254,
        15.196905,
        10.656576,
        21.494133,
        22.761848,
        19.585612,
        22.843220,
        15.911071,
        20.953345,
        22.681723,
        19.391355,
        15.748333,
        10.002256,
        21.294325,
        22.634677,
        19.515170,
        17.858554,
        9.908810,
    ]
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_B"]
            == pytest.approx(expected_PAE_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )
