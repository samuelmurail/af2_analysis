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

    assert my_data.format == "colabfold_1.5"
    assert len(my_data.df) == 40
    print(my_data.df.columns)
    assert (
        my_data.df.columns
        == np.array(
            [
                "query",
                "seed",
                "model",
                "weight",
                "recycle",
                "pLDDT",
                "pTM",
                "ipTM",
                "ranking_confidence",
                "pdb",
                "relaxed_pdb",
                "json",
            ]
        )
    ).all()

    query = my_data.df.iloc[0]["query"]

    assert my_data.chain_length[query] == [42, 42]
    assert my_data.chains[query] == ["A", "B"]

    # There should be only 5 relaxed structures

    relaxed_num = sum(my_data.df["relaxed_pdb"].notna())

    assert relaxed_num == 5
    assert list(my_data.df["recycle"]) == [
        9,
        16,
        3,
        5,
        5,
        14,
        15,
        15,
        48,
        12,
        7,
        14,
        7,
        8,
        9,
        4,
        18,
        6,
        10,
        8,
        11,
        14,
        7,
        11,
        9,
        11,
        7,
        7,
        19,
        9,
        30,
        7,
        5,
        9,
        7,
        12,
        7,
        6,
        6,
        6,
    ]

    assert list(my_data.df["ipTM"]) == [
        0.0812,
        0.0685,
        0.161,
        0.158,
        0.541,
        0.117,
        0.0698,
        0.239,
        0.0648,
        0.331,
        0.0789,
        0.0815,
        0.145,
        0.306,
        0.604,
        0.0997,
        0.0662,
        0.143,
        0.219,
        0.589,
        0.0794,
        0.0684,
        0.15,
        0.299,
        0.559,
        0.0797,
        0.0662,
        0.147,
        0.0609,
        0.318,
        0.0964,
        0.0683,
        0.151,
        0.274,
        0.584,
        0.0776,
        0.0693,
        0.14,
        0.199,
        0.598,
    ]

def test_af3_webserver():
    data_path = os.path.join(TEST_FILE_PATH, "fold_2024_07_01_12_14_prot_dna_zn")

    my_data = af2_analysis.Data(data_path)
    analysis.pdockq(my_data)

    expected_pdockq = [
        0.1984, 0.2221, 0.2478, 0.3059, 0.2736
    ]
    assert 2.2 == pytest.approx(2.3, 0.1)
    print([round(i, 4) for i in my_data.df["pdockq"]])
    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]["pdockq"] == pytest.approx(expected_pdockq[i], precision)
            for i in range(len(my_data.df))
        ]
    )
    # TO FIX !!!

    #analysis.pdockq2(my_data)

    analysis.LIS_matrix(my_data)
    expected_LIS_0 = [[0.70403, 0.59714, 0.62242, 0.67410, 0.71366, 0.69927],
                      [0.20782, 0.65753, 0.46597, 0.06869, 0.19638, 0.11468],
                      [0.32243, 0.54838, 0.70598, 0.20183, 0.30740, 0.38108],
                      [0.59661, 0.59681, 0.57492, 0.93666, 0.57333, 0.43416],
                      [0.65602, 0.64151, 0.64242, 0.58833, 0.93666, 0.55083],
                      [0.61730, 0.54878, 0.63575, 0.3825 , 0.52916, 0.93666],]
    print(np.array(my_data.df["LIS"][0]))
    np.testing.assert_allclose(
        np.array(my_data.df["LIS"][0]),
        np.array(expected_LIS_0),
        atol=precision)

    analysis.inter_chain_pae(my_data)

    expected_PAE_A_B = [5.0452, 3.9899, 4.3717, 3.928, 3.5882]
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_B"]
            == pytest.approx(expected_PAE_A_B[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    expected_PAE_A_E = [3.9401, 3.1581, 3.5433, 3.5354, 3.1004]
    print([round(i, 4) for i in my_data.df["PAE_A_E"]])
    assert np.all(
        [
            my_data.df.iloc[i]["PAE_A_E"]
            == pytest.approx(expected_PAE_A_E[i], precision)
            for i in range(len(my_data.df))
        ]
    )

