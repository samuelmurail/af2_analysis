import numpy as np

# Autorship information
__author__ = "Alaa Reguei"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.2"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"

def compute_LIS_matrix(
        coor,
        pae_array,
        pae_cutoff=12.0,):
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
    chain_lens = [coor.select_atoms(f"name CA and chain {chain}").len for chain in chain_list]
    chain_len_sums = np.cumsum([0] + chain_lens)

    LIS_array = np.zeros((len(chain_lens),len(chain_lens)))

    trans_matrix = np.zeros_like(pae_array)
    mask = pae_array < pae_cutoff
    trans_matrix[mask] = 1 - pae_array[mask]/pae_cutoff
    
    for i in range(len(chain_list)):

        i_start = chain_len_sums[i]
        i_end = chain_len_sums[i+1]
        
        for j in range(len(chain_list)):
            
            j_start = chain_len_sums[j]
            j_end = chain_len_sums[j+1]

            submatrix = trans_matrix[i_start: i_end, j_start:j_end]

            mean_lis = submatrix[submatrix > 0].mean()
            if not np.isnan(mean_lis):
                LIS_array[i, j] = mean_lis

    return(LIS_array)