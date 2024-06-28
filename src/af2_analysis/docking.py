import numpy as np
import pdb_numpy

from tqdm.auto import tqdm
from . import data

"""
The module contains functions to extract and compute docking scores.

..Warning:
    The ligand chain is assumed to be the last chain in the list of chains.
"""


def extract_pae_pep(my_data, fun=np.mean, verbose=True):
    """Extract the PAE score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the PAE scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    my_data.extract_inter_chain_pae(verbose=verbose)

    pep_rec_pae_list = []
    rec_pep_pae_list = []

    disable = False if verbose else True

    for i, (query, json) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["json"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        pae = data.get_pae(json)

        # print(f"{cum_sum_chain}  {cum_sum_chain[-2]}:{cum_sum_chain[-1]}")
        rec_pep_pae = fun(
            pae[0 : cum_sum_chain[-2], cum_sum_chain[-2] : cum_sum_chain[-1]]
        )
        pep_rec_pae = fun(
            pae[cum_sum_chain[-2] : cum_sum_chain[-1], 0 : cum_sum_chain[-2]]
        )

        pep_rec_pae_list.append(pep_rec_pae)
        rec_pep_pae_list.append(rec_pep_pae)

    my_data.df.loc[:, "PAE_pep_rec"] = pep_rec_pae_list
    my_data.df.loc[:, "PAE_rec_pep"] = rec_pep_pae_list


def extract_plddt_pep(my_data, fun=np.mean, verbose=True):
    """Extract the pLDDT score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the pLDDT scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """
    pep_plddt_list = []

    disable = False if verbose else True

    for i, (query, pdb) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["pdb"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        plddt = my_data.get_plddt(i)
        pep_plddt_list.append(fun(plddt[cum_sum_chain[-2] : cum_sum_chain[-1]]))

    my_data.df.loc[:, "plddt_pep"] = pep_plddt_list


def extract_plddt_contact_pep(my_data, fun=np.mean, cutoff=8.0, verbose=True):
    """Extract the pLDDT score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    fun : function
        function to apply to the pLDDT scores

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """
    lig_plddt_list = []
    rec_plddt_list = []

    disable = False if verbose else True

    for i, (query, pdb) in tqdm(
        enumerate(zip(my_data.df["query"], my_data.df["pdb"])),
        total=len(my_data.df),
        disable=disable,
    ):
        chain_length = my_data.chain_length[query]
        chains = my_data.chains[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        model = pdb_numpy.Coor(pdb)
        model_CA = model.select_atoms("name CA")
        contact_lig = model_CA.select_atoms(
            f"chain {chains[-1]} and within {cutoff} of chain {' '.join(chains[:-1])}"
        )
        contact_rec = model_CA.select_atoms(
            f"chain {' '.join(chains[:-1])} and within {cutoff} of chain {chains[-1]}"
        )

        if contact_lig.len > 0:
            lig_plddt_list.append(fun(contact_lig.beta))
        else:
            lig_plddt_list.append(0)

        if contact_rec.len > 0:
            rec_plddt_list.append(fun(contact_rec.beta))
        else:
            rec_plddt_list.append(0)

    my_data.df.loc[:, "plddt_contact_lig"] = lig_plddt_list
    my_data.df.loc[:, "plddt_contact_rec"] = rec_plddt_list


def compute_LIS_pep(my_data, pae_cutoff=12.0, fun=np.max, verbose=True):
    """Compute the LIS score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A
    fun : function
        function to apply to the LIS matrix

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.

    """

    my_data.compute_LIS_matrix(pae_cutoff=pae_cutoff, verbose=verbose)

    pep_LIS_list = []
    pep_LIS2_list = []

    for query, LIS in zip(my_data.df["query"], my_data.df["LIS"]):
        chain_num = len(my_data.chains[query])

        pep_LIS_list.append(fun(LIS[0 : chain_num - 1, chain_num - 1]))
        pep_LIS2_list.append(fun(LIS[chain_num - 1, 0 : chain_num - 1]))

    my_data.df.loc[:, "LIS_rec_pep"] = pep_LIS2_list
    my_data.df.loc[:, "LIS_pep_rec"] = pep_LIS_list


def compute_pdockq2_lig(my_data, verbose=True):
    """Compute the LIS score for the peptide-peptide interface.

    Parameters
    ----------
    my_data : AF2Data
        object containing the data
    pae_cutoff : float
        cutoff for native contacts, default is 8.0 A
    fun : function
        function to apply to the LIS matrix

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.

    """

    my_data.compute_pdockq2(verbose=verbose)

    old_query = ""
    pdockq2_list = []

    for index, row in my_data.df.iterrows():
        if row["query"] != old_query:
            old_query = row["query"]
            chains = my_data.chains[old_query]
            lig_chain = chains[-1]
            rec_chains = chains[:-1]

        pdockq2_list.append(row[f"pdockq2_{lig_chain}"])

    my_data.df.loc[:, "pdockq2_lig"] = pdockq2_list
