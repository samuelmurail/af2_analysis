#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest

import af2_analysis
from af2_analysis import analysis
from .data_files import TEST_FILE_PATH

def test_cf_1_5_5_write_read_csv(tmp_path):
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")
    my_data = af2_analysis.Data(data_path)
    my_data.export_csv(os.path.join(tmp_path, "test.csv"))

    my_data2 = af2_analysis.Data(csv=os.path.join(tmp_path, "test.csv"))

    assert my_data2.format == "csv"
    assert len(my_data2.df) == 40
    print(my_data2.df.columns)
    assert (
        my_data2.df.columns
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

    query = my_data2.df.iloc[0]["query"]

    assert my_data2.chain_length[query] == [42, 42]
    assert my_data2.chains[query] == ["A", "B"]

    # There should be only 5 relaxed structures

    relaxed_num = sum(my_data2.df["relaxed_pdb"].notna())

    assert relaxed_num == 5
    assert list(my_data2.df["recycle"]) == [
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

    assert list(my_data2.df["ipTM"]) == [
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
