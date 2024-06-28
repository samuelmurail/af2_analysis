#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pytest

import af2_analysis
from af2_analysis import analysis
from .data_files import TEST_FILE_PATH

def test_cf_1_5_5_relax():

    data_path = os.path.join(TEST_FILE_PATH, 'beta_amyloid_dimer_cf_1.5.5')

    my_data = af2_analysis.Data(data_path)
    analysis.pdockq(my_data)

    expected_pdockq = [
        0.03322, 0.02051, 0.03, 0.02971, 0.12094,
        0.02849, 0.02246, 0.04845, 0.0198, 0.07149,
        0.03319, 0.02379, 0.02763, 0.05582, 0.13833,
        0.02418, 0.02112, 0.02521, 0.04189, 0.14153,
        0.02952, 0.02123, 0.02851, 0.05237, 0.13703,
        0.02909, 0.02042, 0.02843, 0.02069, 0.08227,
        0.02312, 0.02025, 0.02824, 0.05085, 0.13916,
        0.0294, 0.02063, 0.02542, 0.03623, 0.14257]
    assert 2.2 == pytest.approx(2.3, 0.1)

    precision = 0.001
    assert np.all([my_data.df.iloc[i]['pdockq'] == pytest.approx(expected_pdockq[i], precision) for i in range(len(my_data.df))])

    analysis.pdockq2(my_data)
    print([round(i,5) for i in my_data.df['pdockq2_B']])

    expected_pdockq2 = [0.01199, 0.01062, 0.01429, 0.01636, 0.16844, 0.01661, 0.0115, 0.02388, 0.00895, 0.02964, 0.01174, 0.01256, 0.01341, 0.05058, 0.20396, 0.01398, 0.01116, 0.01292, 0.02662, 0.19487, 0.01218, 0.01122, 0.01373, 0.05064, 0.18241, 0.01205, 0.01135, 0.01352, 0.01025, 0.02231, 0.01376, 0.0106, 0.01412, 0.04422, 0.19023, 0.01247, 0.01078, 0.01258, 0.02291, 0.20544]
    assert np.all([my_data.df.iloc[i]['pdockq2_B'] == pytest.approx(expected_pdockq2[i], precision) for i in range(len(my_data.df))])

    analysis.LIS_matrix(my_data)
    analysis.inter_chain_pae(my_data)  

