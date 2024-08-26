#!/usr/bin/env python3
# coding: utf-8

import os
import pytest
import numpy as np
import af2_analysis
from af2_analysis import clustering
from .data_files import TEST_FILE_PATH
from unittest.mock import patch 

@patch("matplotlib.pyplot.show")
def test_cf_1_5_5_relax(mock_show):
    data_path = os.path.join(TEST_FILE_PATH, "beta_amyloid_dimer_cf_1.5.5")

    my_data = af2_analysis.Data(data_path)

    assert 'cluster' not in my_data.df.columns
    clustering.hierarchical(my_data.df, threshold=0.3)
    assert 'cluster' in my_data.df.columns

    clust_num = len(my_data.df['cluster'].unique())
    assert clust_num == 10

    assert list(my_data.df['cluster']) == [8, 1, 6, 7, 7, 2, 1, 7, 9, 5, 8, 4, 6, 7, 7, 3, 1, 6, 7, 7, 8, 1, 6, 7, 7, 8, 1, 6, 10, 5, 2, 3, 6, 7, 7, 8, 1, 6, 7, 7]

    clustering.hierarchical(my_data.df)
    clust_num = len(my_data.df['cluster'].unique())
    assert clust_num == 12
    print([i for i in my_data.df['cluster']])
    assert list(my_data.df['cluster']) == [10, 2, 7, 8, 8, 3, 1, 8, 11, 6, 10, 5, 7, 8, 8, 4, 1, 7, 8, 8, 9, 1, 7, 8, 8, 9, 2, 7, 12, 6, 3, 4, 7, 8, 8, 9, 1, 7, 8, 8]

    clustering.compute_pc(my_data.df)
    assert 'PC1' in my_data.df.columns
    assert 'PC2' in my_data.df.columns
    assert 'PC3' not in my_data.df.columns
    clustering.compute_pc(my_data.df, n_components=3)
    assert 'PC3' in my_data.df.columns
    print([round(i, 5) for i in my_data.df['PC1']])

    expected_PC1 = [-95.77812, -94.68417, 75.86093, 65.48964, 83.74089, -60.82929, -110.32781, 79.62047, -95.07062, 92.42451, -101.76761, -103.73307, 74.16488, 78.19159, 85.69576, -103.47705, -110.98566, 75.02448, 72.04246, 85.72818, -94.67249, -108.25562, 74.5274, 74.19903, 83.40501, -95.68722, -95.98197, 73.73635, -90.01166, 90.86801, -53.2237, -107.45218, 74.89252, 73.73373, 86.01567, -97.25223, -108.42328, 77.9843, 70.42554, 79.8424]

    precision = 0.001
    assert np.all(
        [
            my_data.df.iloc[i]['PC1'] == pytest.approx(expected_PC1[i], precision)
            for i in range(len(my_data.df))
        ]
    )

    clustering.plot_pc(my_data.df)