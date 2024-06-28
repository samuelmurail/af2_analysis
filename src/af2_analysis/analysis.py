import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import json

import pdb_numpy

# Autorship information
__author__ = "Alaa Reguei"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.2"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"


def get_pae(json_file):
    """Get the PAE matrix from a json file.

    Parameters
    ----------
    json_file : str
        Path to the json file.

    Returns
    -------
    np.array
        PAE matrix.
    """

    if json_file is None:
        return (None, None)

    with open(json_file) as f:
        local_json = json.load(f)

    if "pae" in local_json:
        pae_array = np.array(local_json["pae"])
    elif "predicted_aligned_error" in local_json[0]:
        pae_array = np.array(local_json[0]["predicted_aligned_error"])
    else:
        raise ValueError("No PAE found in the json file.")

    return pae_array


def pdockq(data, verbose=True):
    """
    Compute the pdockq from the pdb file.

    Parameters
    ----------
    data : AF2Data
        object containing the data
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    from pdb_numpy.analysis import compute_pdockQ

    pdockq_list = []

    disable = False if verbose else True

    for pdb in tqdm(data.df["pdb"], total=len(data.df["pdb"]), disable=disable):
        if pdb:
            model = pdb_numpy.Coor(pdb)
            pdockq_list += compute_pdockQ(model)
        else:
            pdockq_list.append(None)

    data.df["pdockq"] = pdockq_list


def mpdockq(data, verbose=True):
    """
    Compute mpdockq from the pdb file.

    """

    from pdb_numpy.analysis import compute_pdockQ

    pdockq_list = []
    disable = False if verbose else True

    for pdb in tqdm(data.df["pdb"], total=len(data.df["pdb"]), disable=disable):
        if pdb is not None and pdb is not np.nan:
            model = pdb_numpy.Coor(pdb)
            pdockq_list += compute_pdockQ(
                model, cutoff=8.0, L=0.728, x0=309.375, k=0.098, b=0.262
            )
        else:
            pdockq_list.append(None)

    data.df.loc[:, "mpdockq"] = pdockq_list


def pdockq2(data, verbose=True):
    r"""
    Compute pdockq2 from the pdb file [1]_.

    .. math::
        pDockQ_2 = \frac{L}{1 + exp [-k*(X_i-X_0)]} + b

    with

    .. math::
        X_i = \langle \frac{1}{1+(\frac{PAE_{int}}{d_0})^2} \rangle - \langle pLDDT \rangle_{int}

    References:

    .. [1] : https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
    """

    from pdb_numpy.analysis import compute_pdockQ2

    pdockq_list = []

    max_chain_num = 0
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num

    for i in range(max_chain_num):
        pdockq_list.append([])

    disable = False if verbose else True

    for pdb, json_path in tqdm(
        zip(data.df["pdb"], data.df["json"]), total=len(data.df["pdb"]), disable=disable
    ):
        if (
            pdb is not None
            and pdb is not np.nan
            and json_path is not None
            and json_path is not np.nan
        ):
            model = pdb_numpy.Coor(pdb)
            # with open(json_path) as f:
            #     local_json = json.load(f)
            # pae_array = np.array(local_json["pae"])
            pae_array = get_pae(json_path)

            pdockq2 = compute_pdockQ2(model, pae_array)

            for i in range(max_chain_num):
                if i < len(pdockq2):
                    pdockq_list[i].append(pdockq2[i][0])
                else:
                    pdockq_list[i].append(None)

        else:
            for i in range(max_chain_num):
                pdockq_list[i].append(None)

    # print(pdockq_list)
    for i in range(max_chain_num):
        data.df.loc[:, f"pdockq2_{chr(65+i)}"] = pdockq_list[i]


def inter_chain_pae(data, fun=np.mean, verbose=True):
    """Read the PAE matrix and extract the average inter chain PAE.

    Parameters
    ----------
    data : AF2Data
        object containing the data
    fun : function
        function to apply to the PAE scores
    verbose : bool
        print progress bar

    Returns
    -------
    None
    """
    pae_list = []

    disable = False if verbose else True

    for query, json_path in tqdm(
        zip(data.df["query"], data.df["json"]),
        total=len(data.df["json"]),
        disable=disable,
    ):
        if json_path is not None and json_path is not np.nan:
            with open(json_path) as f:
                local_json = json.load(f)
            pae_array = np.array(local_json["pae"])

            chain_lens = data.chain_length[query]
            chain_len_sums = np.cumsum([0] + chain_lens)
            # pae_chain_array = np.empty((len(chain_lens), len(chain_lens)))
            chain_ids = data.chains[query]

            pae_dict = {}

            for i in range(len(chain_lens)):
                for j in range(len(chain_lens)):
                    pae_val = fun(
                        pae_array[
                            chain_len_sums[i] : chain_len_sums[i + 1],
                            chain_len_sums[j] : chain_len_sums[j + 1],
                        ]
                    )
                    pae_dict[f"PAE_{chain_ids[i]}_{chain_ids[j]}"] = pae_val

            pae_list.append(pae_dict)
        else:
            pae_list.append({})

    pae_df = pd.DataFrame(pae_list)

    for col in pae_df.columns:
        data.df.loc[:, col] = pae_df.loc[:, col].to_numpy()


def compute_LIS_matrix(
    coor,
    pae_array,
    pae_cutoff=12.0,
):
    r"""Compute the LIS score as define in [1]_.

    Implementation was inspired from implementation in https://github.com/flyark/AFM-LIS

    Parameters
    ----------
    coor : Coor
        object containing the coordinates of the model
    pae_array : np.array
        array of predicted PAE
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A

    Returns
    -------
    list
        LIS scores

    References
    ----------

    .. [1] Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. "Enhanced
        Protein-Protein Interaction Discovery via AlphaFold-Multimer" bioRxiv (2024).
        https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1

    """

    chain_list = np.unique(coor.chain)
    chain_lens = [
        coor.select_atoms(f"name CA and chain {chain}").len for chain in chain_list
    ]
    chain_len_sums = np.cumsum([0] + chain_lens)

    LIS_array = np.zeros((len(chain_lens), len(chain_lens)))

    trans_matrix = np.zeros_like(pae_array)
    mask = pae_array < pae_cutoff
    trans_matrix[mask] = 1 - pae_array[mask] / pae_cutoff

    for i in range(len(chain_list)):
        i_start = chain_len_sums[i]
        i_end = chain_len_sums[i + 1]

        for j in range(len(chain_list)):
            j_start = chain_len_sums[j]
            j_end = chain_len_sums[j + 1]

            submatrix = trans_matrix[i_start:i_end, j_start:j_end]

            if np.any(submatrix > 0):
                LIS_array[i, j] = submatrix[submatrix > 0].mean()

    return LIS_array


def LIS_matrix(data, pae_cutoff=12.0, verbose=True):
    """
    Compute the LIS score as define in [2]_.

    Implementation was inspired from implementation in:

    .. [2] https://github.com/flyark/AFM-LIS

    Parameters
    ----------
    data : AF2Data
        object containing the data
    pae_cutoff : float
        cutoff for PAE matrix values, default is 12.0 A
    verbose : bool
        print progress bar

    Returns
    -------
    None
        The dataframe is modified in place.
    """
    LIS_matrix_list = []

    max_chain_num = 0
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num

    disable = False if verbose else True

    for pdb, json_path in tqdm(
        zip(data.df["pdb"], data.df["json"]), total=len(data.df["pdb"]), disable=disable
    ):
        if (
            pdb is not None
            and pdb is not np.nan
            and json_path is not None
            and json_path is not np.nan
        ):
            model = pdb_numpy.Coor(pdb)
            pae_array = get_pae(json_path)
            LIS_matrix = compute_LIS_matrix(model, pae_array, pae_cutoff)
            LIS_matrix_list.append(LIS_matrix)
        else:
            LIS_matrix_list.append(None)

    data.df.loc[:, "LIS"] = LIS_matrix_list


def piTM(data, verbose=True):
    r"""Compute the piTM score as define in [3]_.

    .. math::
        piTM = \max_{i \in \mathcal{I}} \frac{1}{I} \sum_{j \in \mathcal{I}}  \frac{1}{1 + [\langle e_{ij} \rangle / d_0 (I)]^2}

    with:

    .. math::
        d_0(I) = \begin{cases} 1.25 \sqrt[3]{I -15} -1.8\text{,} & \text{if } I \geq 22 \\ 0.02 I \text{,} & \text{if } I < 22  \end{cases}


    Implementation was inspired from `predicted_tm_score_v1()` in https://github.com/FreshAirTonight/af2complex/blob/main/src/alphafold/common/confidence.py

    References
    ----------
    .. [3] Mu Gao, Davi Nakajima An, Jerry M. Parks & Jeffrey Skolnick.
        AF2Complex predicts direct physical interactions in multimeric proteins with deep learning
        *Nature Communications*. volume 13, Article number: 1744 (2022).
        https://www.nature.com/articles/s41467-022-29394-2

    .. warning::
        IT IS NOT WORKING !!
    """

    from pdb_numpy.analysis import compute_piTM

    piTM_chain_list = []
    piTM_list = []

    max_chain_num = 0
    for query in data.chains:
        chain_num = len(data.chains[query])
        if chain_num > max_chain_num:
            max_chain_num = chain_num

    for i in range(max_chain_num):
        piTM_chain_list.append([])

    disable = False if verbose else True

    for pdb, json_path in tqdm(
        zip(data.df["pdb"], data.df["json"]), total=len(data.df["pdb"]), disable=disable
    ):
        # print(pdb, json_path)
        if (
            pdb is not None
            and pdb is not np.nan
            and json_path is not None
            and json_path is not np.nan
        ):
            model = pdb_numpy.Coor(pdb)
            with open(json_path) as f:
                local_json = json.load(f)
            pae_array = np.array(local_json["pae"])

            piTM, piTM_chain = compute_piTM(model, pae_array)
            # print(piTM, piTM_chain)

            piTM_list.append(piTM[0])

            for i in range(max_chain_num):
                # print(i, len(piTM_chain))
                if i < len(piTM_chain):
                    piTM_chain_list[i].append(piTM_chain[i][0])
                else:
                    piTM_chain_list[i].append(None)

        else:
            piTM_list.append(None)
            for i in range(max_chain_num):
                piTM_chain_list[i].append(None)

    # print(piTM_chain_list)
    data.df.loc[:, "piTM"] = piTM_list
    for i in range(max_chain_num):
        data.df.loc[:, f"piTM_{chr(65+i)}"] = piTM_chain_list[i]
