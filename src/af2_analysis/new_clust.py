import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import sklearn.manifold
import sklearn.decomposition
from MDAnalysis.analysis import align, diffusionmap, pca
from MDAnalysis.coordinates.memory import MemoryReader
import pandas as pd
import pdb_numpy
import os
import shutil
import logging


# Autorship information
__author__ = "Alaa Reguei, Samuel Murail"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.2"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"


"""
The module contains functions to cluster AlphaFold models.
"""

# Logging
logger = logging.getLogger(__name__)


def read_numerous_pdb(pdb_files, batch_size=1000):
    """Read a large number of PDB files in batches and combine
    them into a single MDAnalysis Universe.

    Discussed in:
    https://github.com/MDAnalysis/mdanalysis/issues/4590

    Parameters
    ----------
    pdb_files : list of str
        List of file paths to the PDB files to be read.
    batch_size : int, optional
        Number of PDB files to read in each batch. Default is 1000.

    Returns
    -------
    MDAnalysis.Universe
        A single MDAnalysis Universe containing the combined frames from all the PDB files.

    Notes
    -----
    - This function reads PDB files in batches to avoid memory issues.
    - Each batch of PDB files is loaded into a temporary Universe, and the positions of each frame
      are stored in a list.
    - The list of frames is then converted into a numpy array and used to create a new Universe with
      a MemoryReader, combining all the frames.

    Example
    -------
    >>> pdb_files = ['file1.pdb', 'file2.pdb', ...]
    >>> combined_universe = read_numerous_pdb(pdb_files, batch_size=1000)
    >>> print(combined_universe)
    """

    all_frames = []

    if pdb_files[0].endswith(".cif"):
        model = pdb_numpy.Coor(pdb_files[0])
        for file in pdb_files[1:]:
            local_model = pdb_numpy.Coor(file)
            model.models.append(local_model.models[0])
        model.write("tmp.pdb", overwrite=True)
        return mda.Universe("tmp.pdb", "tmp.pdb")

    for i in range(0, len(pdb_files), batch_size):
        # print(f"Reading frames {i:5} to {i+batch_size:5}, total : {len(pdb_files[i:i+batch_size])} frames")
        local_u = mda.Universe(pdb_files[0], pdb_files[i : i + batch_size])
        for ts in local_u.trajectory:
            all_frames.append(ts.positions.copy())
        del local_u

    # print("Convert to numpy")
    frames_array = np.array(all_frames)
    del all_frames

    # print(frames_array.shape)
    return mda.Universe(pdb_files[0], frames_array, format=MemoryReader, order="fac")

def scale(rms, d0=8.5):
    """Scaling RMS values.
    To avoid the problem of arbitrarily large RMS values that are essentially equally bad,
    this function scales RMS values using the inverse square scaling technique adapted
    from the S-score formula.

    [Sankar Basu, Björn Wallner https://doi.org/10.1371/journal.pone.0161879]

    Parameters
    ----------
    rms : float or numpy.ndarray
        RMS values as single values or matrix.

    d0 : float
        Scaling factor, fixed to 8.5 for LRMS.

    Returns
    -------
    rms_scale : float or numpy.ndarray
        Scaled RMS values as single values or matrix.
    """

    rms_scale = 1 / (1 + (rms / d0) ** 2)
    return rms_scale

def compute_distance_matrix(
    pdb_files,
    align_selection=None,
    distance_selection=None,
):
    """Compute distance matrix.
    
    This function computes the distance matrix of the different peptide chains.
    
    Parameters
    ----------
    pdb_files : list of str
        List of file paths to the PDB files to be read.
    align_selection : str, optional
        The selection string to align the protein chains (default is None or "backbone").
    distance_selection : str, optional
        The selection string to compute the distance matrix (default is None or "backbone").
    
    Returns
    -------
    None
        The function modifies the input DataFrame by appending the DistanceMatrix column, but does not return a value.
    """

    if align_selection is None and distance_selection is None:
        distance_selection = "backbone"
        align_selection = "backbone"
    elif align_selection is None:
        align_selection = distance_selection
    elif distance_selection is None:
        distance_selection = align_selection
    
    u = read_numerous_pdb(pdb_files)
    assert u.trajectory.n_frames == len(
        pdb_files
    ), f"Number of frames {u.trajectory.n_frames} different from number of files {len(files)}"

    logger.info("align structures")
    align.AlignTraj(
        u, u, select=align_selection, in_memory=True
    ).run(verbose=True)

    logger.info("Compute distance Matrix")
    matrix = diffusionmap.DistanceMatrix(
        u,
        select=distance_selection,
    ).run(verbose=True)

    return matrix.results.dist_matrix

    
def hierarchical(
    df,
    threshold=0.2,
    align_selection=None,
    distance_selection=None,
    show_dendrogram=True,
    rmsd_scale=False,
):
    """Clustering of AlphaFold Protein-Peptide Complex results.

    After checking for the absence of missing values, the function starts by aligning the protein
    chains of different predicted models of each PDB before characterizing the different
    protein-peptide contact residues. The search for contact residues is done within 4 angstroms
    by default. These residues will be used to compute the distance matrix of the different
    peptide chains, where the output array will be scaled using the `scale` function before
    calculating the new matrix of 1-scaled matrix. This new matrix will be used to perform
    hierarchical ascending classification and characterize clusters by defining a cutoff threshold.

    Optionally, this function can also plot the dendrogram of each PDB and the clusters distribution
    plot using the `clusters_distribution`.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing AlphaFold data for the protein-peptide complex.
    threshold : float
        The value where to cut the dendrogram in order to define clusters.
    align_selection : str or dict, optional
        The selection string to align the protein chains (default is None).
    distance_selection : str or dict, optional
        The selection string to compute the distance matrix (default is None).
    show_dendrogram : bool, optional
        Whether to plot the different dendrograms of each PDB (default is `True`).
    rmsd_scale : bool, optional
        Whether to scale the RMS values using the `scale` function (default is `False`).

    Returns
    -------
    None
        The function modifies the input DataFrame by appending the Cluster column, but does not return a value.
    """

    cluster_list = []
    MDS_1 = []
    MDS_2 = []

    query_list = df["query"].unique().tolist()

    if align_selection is None:
        align_selection = {pdb: "backbone" for pdb in query_list}
    elif isinstance(align_selection, str):
        align_selection = {pdb: align_selection for pdb in query_list}

    if distance_selection is None:
        distance_selection = {pdb: "backbone" for pdb in query_list}
    elif isinstance(distance_selection, str):
        distance_selection = {pdb: distance_selection for pdb in query_list}

    for pdb in query_list:
        files = [file for file in df[df["query"] == pdb]["pdb"] if not pd.isnull(file)]

        # Check that only the last df rows are missing
        null_number = sum(pd.isnull(df[df["query"] == pdb]["pdb"]))
        assert null_number == sum(
            pd.isnull(df[df["query"] == pdb]["pdb"][-null_number:])
        ), f"Missing pdb data in the middle of the query {pdb}"

        logger.info("Read all structures")
        dist_matrix = compute_distance_matrix(
            files,
            align_selection=align_selection[pdb],
            distance_selection=distance_selection[pdb],)

        logger.info(f"Max RMSD is {np.max(dist_matrix):.2f} A")
        if rmsd_scale:
            dist_matrix = 1 - scale(dist_matrix)

        logger.info("Compute Linkage clustering")
        h, _ = dist_matrix.shape
        Z = linkage(dist_matrix[np.triu_indices(h, 1)], method="average")
        cluster_list += fcluster(Z, float(threshold), criterion="distance").tolist()

        logger.info(f"{len(np.unique(cluster_list))} clusters founded for {pdb}")

        # Add Multidimensional scaling coordinates from the distance matrix
        mds = sklearn.manifold.MDS(dissimilarity='precomputed', n_components=2)
        coordinates = mds.fit_transform(dist_matrix)

        # PCA analysis don't add any information:
        #pca = sklearn.decomposition.PCA(n_components=2)  # n_components defines how many principal components you want
        #pca_result = pca.fit_transform(coordinates)
        #coordinates = pca_result
        #print(pca_result)

        if show_dendrogram:
            # plot the dendrogram with the threshold line
            plt.figure()
            dendrogram(Z, color_threshold=float(threshold))
            plt.axhline(float(threshold), color="k", ls="--")
            plt.title("Hierarchical Cluster Dendrogram -{}".format(pdb))
            plt.xlabel("Data Point Indexes")
            plt.ylabel("Distance")
            plt.show()

    df["cluster"] = cluster_list + null_number * [None]
    df["cluster"] = df["cluster"].astype('category')
    df["MDS 1"] = list(coordinates.T[0]) + null_number * [None]
    df["MDS 2"] = list(coordinates.T[1]) + null_number * [None]


def compute_pc(df, n_components=2):
    """Compute Principal Components for Alphafold Protein-Peptide Complex.

    This function computes the Principal Components (PCs) for a protein-peptide complex predicted by
    Alphafold. After checking for the absence of missing values, the function aligns the protein
    chains before calculating the PCs using the PCA module from the MDAnalysis package. The results
    are then updated in the input DataFrame. By default, this function computes the first two
    principal components.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing Alphafold data for the protein-peptide complex.
    n_components : int, optional
        Number of principal components to compute (default is 2).

    Returns
    -------
    None
        The function modifies the input DataFrame by appending the PC columns, but does not return a value.
    """

    pc_list = []
    for i in range(n_components):
        pc_list.append([])

    query_list = df["query"].unique().tolist()

    for pdb in query_list:
        files = [file for file in df[df["query"] == pdb]["pdb"] if not pd.isnull(file)]
        # Check that only the last df row are missing
        null_number = sum(pd.isnull(df[df["query"] == pdb]["pdb"]))
        assert null_number == sum(
            pd.isnull(df[df["query"] == pdb]["pdb"][-null_number:])
        ), f"Missing pdb data in the middle of the query {pdb}"
        u = read_numerous_pdb(files)
        chain_pep_value = chain_pep(files[0])

        align.AlignTraj(
            u, u, select=f"backbone and not chainID {chain_pep_value}", in_memory=True
        ).run()
        pc = pca.PCA(
            u,
            select=(f"backbone and chainID {chain_pep_value}"),
            align=True,
            mean=None,
            n_components=n_components,
        ).run()
        backbone = u.select_atoms(f"backbone and chainID {chain_pep_value}")
        transformed = pc.transform(backbone, n_components)
        for i in range(n_components):
            pc_list[i].extend(transformed[:, i].tolist())

    for i in range(n_components):
        df[f"PC{i+1}"] = pc_list[i] + null_number * [np.nan]
