import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from MDAnalysis.analysis import align, diffusionmap, pca
from MDAnalysis.coordinates.memory import MemoryReader
import pandas as pd
import pdb_numpy
import os
import shutil


# Autorship information
__author__ = "Alaa Reguei"
__copyright__ = "Copyright 2023, RPBS"
__credits__ = ["Samuel Murail", "Alaa Reguei"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.2"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"


def chain_pep(pdb_file):
    """Search for peptide chain ID.

    This function searches for the peptide chain ID of a protein-peptide complex by using
    the PDB numpy package.

    For each PDB file, the local coordinates will be loaded, extracting the chains of the
    complex and limiting to the last one.

    Parameters
    ----------
    pdb_file : str
        Path to the protein-peptide complex PDB file.

    Returns
    -------
    chain : str
        The peptide chain ID.
    """
    local_coor = pdb_numpy.Coor(pdb_file)
    chains = np.unique(local_coor.chain)

    # Get the smallest chain
    min_num = len(local_coor.chain)
    for chain in chains:
        if len(local_coor.chain[local_coor.chain == chain]) < min_num:
            min_num = len(local_coor.chain[local_coor.chain == chain])
            chain_value = chain

    return chain


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


def hierarchical(
    df,
    threshold=0.2,
    contact_cutoff=4.0,
    show_dendrogram=True,
    show_cluster_distribution=False,
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
    show_dendrogram : bool, optional
        Whether to plot the different dendrograms of each PDB (default is `True`).
    show_cluster_distribution : bool, optional
        Whether to plot the clusters distribution of all queries (default is `False`).

    Returns
    -------
    None
        The function modifies the input DataFrame by appending the Cluster column, but does not return a value.
    """

    clust_threshold = []
    query_list = df["query"].unique().tolist()

    for pdb in query_list:
        files = [file for file in df[df["query"] == pdb]["pdb"] if not pd.isnull(file)]

        # Check that only the last df rows are missing
        null_number = sum(pd.isnull(df[df["query"] == pdb]["pdb"]))
        assert null_number == sum(
            pd.isnull(df[df["query"] == pdb]["pdb"][-null_number:])
        ), f"Missing pdb data in the middle of the query {pdb}"

        print("Read all structures")
        u = read_numerous_pdb(files)
        assert u.trajectory.n_frames == len(
            files
        ), f"Number of frames {u.trajectory.n_frames} different from number of files {len(files)}"
        chain_pep_value = chain_pep(files[0])
        print(f"Peptide chain is :{chain_pep_value}")

        print("align structures")
        align.AlignTraj(
            u, u, select=f"backbone and not chainID {chain_pep_value}", in_memory=True
        ).run(verbose=True)

        print("Extract contact residues")
        cutoff = contact_cutoff
        peptide_contact = u.select_atoms(
            f"chainID {chain_pep_value} and around {cutoff} not chainID {chain_pep_value}"
        )
        resid_contact_list = []
        for _ in u.trajectory:
            for residue in peptide_contact.groupby("residues"):
                resid_contact_list.append(residue.resnum)

        resid_contact_list = set(resid_contact_list)
        print(f"Contact residues : {resid_contact_list}")

        print("Compute distance Matrix")
        # matrix = diffusionmap.DistanceMatrix(u, select=f"backbone and chainID {chain_pep_value}").run(verbose=True)
        matrix = diffusionmap.DistanceMatrix(
            u,
            select=f'chainID {chain_pep_value} and resnum {" ".join([str(res) for res in resid_contact_list])} and backbone',
        ).run(verbose=True)

        print(f"Max RMSD is {np.max(matrix.results.dist_matrix):.2f} A")
        dist = 1 - scale(matrix.results.dist_matrix)
        h, w = dist.shape

        print("Compute Linkage clustering")
        Z = linkage(dist[np.triu_indices(h, 1)], method="average")
        clust_threshold.extend(
            fcluster(Z, float(threshold), criterion="distance").tolist()
        )

        print(f"{len(np.unique(clust_threshold))} clusters founded for {pdb}")

        if show_dendrogram:
            # plot the dendrogram with the threshold line
            plt.figure(figsize=(15, 5))
            dendrogram(Z, color_threshold=float(threshold))
            plt.axhline(float(threshold), color="k", ls="--")
            plt.title("Hierarchical Cluster Dendrogram -{}".format(pdb))
            plt.xlabel("Data Point Indexes")
            plt.ylabel("Distance")
            plt.show()

    df["cluster"] = clust_threshold + null_number * [-1]

    if show_cluster_distribution:
        clusters_distribution(df)


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


def clusters_distribution(df):
    """Plotting the clusters distribution.

    This function plots the distribution of clusters of all PDB queries. It starts
    by determining the number of clusters found for each query before preparing
    a DataFrame presenting all found with their corresponding PDB code.

    A histogram of the clusters distribution will be created from the DataFrame.

    Parameters
    ----------
    df : DataFrame
         DataFrame containing Alphafold data.

    Returns
    -------
    None
    """
    query_list = df["query"].unique().tolist()
    n_cluster = []
    for pdb in query_list:
        sub_df = df[df["query"] == pdb]
        n_cluster.append({"pdb": pdb, "n_clusters": len(np.unique(sub_df["cluster"]))})
    df = pd.DataFrame(n_cluster)
    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=df,
        x="n_clusters",
        palette="Dark2",
        hue="n_clusters",
        discrete=True,
        stat="percent",
    )
    plt.xlabel("Number of Clusters")
    plt.ylabel("Frequency")
    plt.title(" Clusters Distribution")
    plt.xlim(0, 90)
    plt.tight_layout()
    plt.show()


def Cluster_reordering(df):
    """Reordering clusters by size.

    This function proposes a new cluster order based on the size of clusters. Starting from
    the DataFrame cluster column, this function ranks the clusters according to their
    sizes before reordering them in a decreasing way from huge clusters to small ones.

    The new cluster order will be added as a column to this input DataFrame.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing Alphafold data.

    Returns
    -------
    df : DataFrame
        Updated Alphafold DataFrame.
    """

    query_list = df["query"].unique().tolist()
    clust_dict = {}
    convert_dict = {}
    new_cluster = []

    for pdb in query_list:
        sub_df = df[df["query"] == pdb]
        n_clusters = len(np.unique(sub_df["cluster"]))

        for clust in range(1, n_clusters + 1):
            if sum(sub_df["cluster"] == clust) > 0:
                clust_dict[clust] = sum(sub_df["cluster"] == clust)

        order_dict = {
            k: v
            for k, v in sorted(
                clust_dict.items(), key=lambda item: item[1], reverse=True
            )
        }
        for i, key in enumerate(order_dict):
            # print(key)
            if key != -1:
                convert_dict[key] = i + 1

        convert_dict[-1] = -1
        sub_df["New_clusters_order"] = [
            convert_dict[clust] for clust in sub_df["cluster"]
        ]
        # create a list of new cluster to add it to the original dataframe
        new_cluster.extend([convert_dict[clust] for clust in sub_df["cluster"]])

    df["New_clusters_order"] = new_cluster


def get_pdb(df, metric, ascending=False):
    """Sorts and copies top structures for each cluster based on the
    provided metric for each query in the DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing Alphafold data for the protein-peptide complex and the different clusters.
    metric : str
        The metric based on which the sorting will be performed.
    ascending : bool, optional
        Specifies the sorting order. If False (default), structures are sorted in descending order of the metric; if True, in ascending order.

    Returns
    -------
    None

    The function sorts the structures based on the provided metric,
    selecting either the top 5 structures for queries with a single cluster
    or the top 3 structures per cluster for queries with multiple clusters.
    It then copies these top structures into respective folders organized
    by query and cluster.
    """
    Clusters = "clusters"

    if not os.path.exists(Clusters):
        os.makedirs(Clusters)

    for pdb in df["query"].unique().tolist():
        sub_df = df[df["query"] == pdb]
        nclusters = sub_df["cluster"].unique().tolist()
        pdb_folder = os.path.join(Clusters, pdb)

        if not os.path.exists(pdb_folder):
            os.makedirs(pdb_folder)

        if len(nclusters) == 1:
            top_structures = (
                sub_df.sort_values(by=metric, ascending=ascending)
                .iloc[0:5]["pdb"]
                .tolist()
            )
            print(
                f"This structure presents {pdb} {len(nclusters)} clusters according to the chosen threshold "
            )
            # print(f"Top 5 structures for {pdb}: {top_structures}")
            for file in top_structures:
                shutil.copy(file, pdb_folder)

        # If nclusters is greater than 2, take the top 3 structure from each cluster
        else:
            print(
                f"This structure presents {pdb} {len(nclusters)} clusters according to the chosen threshold "
            )
            grouped_df = sub_df.groupby("cluster")

            for i in nclusters:
                group_for_cluster_i = grouped_df.get_group(i)
                top_structure = group_for_cluster_i.sort_values(
                    by=metric, ascending=ascending
                ).iloc[0:3]["pdb"]
                # print(f"Top structure for cluster {i} in {pdb}: {top_structure}")
                cluster_folder = os.path.join(pdb_folder, f"cluster_{i}")

                if not os.path.exists(cluster_folder):
                    os.makedirs(cluster_folder)
                for file in top_structure:
                    shutil.copy(file, cluster_folder)


def compute_pc(df, n_components=2):
    """Compute Principal Components for Alphafold Protein-Peptide Complex.

    This function computes the Principal Components (PCs) for a protein-peptide complex predicted by
    Alphafold. After checking for the absence of missing values, the function aligns the protein
    chains before calculating the PCs using the PCA module from the MDAnalysis package. The results
    are then updated in the input DataFrame. By default, this function computes the first three
    principal components. Optionally, it can also plot the PCs, defaulting to the first two components
    if the `plot_pca` parameter is set to `True`.

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


def plot_pc(df, X="PC1", Y="PC2", show_legend=False, min_clust_num=5, **kwargs):
    """Plotting principal components.

    This function helps in visualizing the pre-computed principal
    components (PCs) and the distribution of clusters for each PDB
    already determined through a clustering function. By default,
    this function plots the first two PCs over a scatterplot. It
    can be directly executed through the `compute_pc` function.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing Alphafold data.
    X : str ,optional
        Name of the column representing the X-axis (default is 'PC1').
    Y : str ,optional
        Name of the column representing the Y-axis (default is 'PC2').
    show_legend : bool, optional
        Whether to show the legend (default is `False`).

    Returns
    -------
    None
    """
    query_list = df["query"].unique().tolist()
    for pdb in query_list:
        sub_df = df[df["query"] == pdb]
        # sub_df["cluster"] = list(sub_df["cluster"])

        # Filter out cluster that are less than min_clust_num
        keep_clust = []
        for cluster in sub_df["cluster"].unique():
            if len(sub_df[sub_df["cluster"] == cluster]) >= min_clust_num:
                keep_clust.append(cluster)

        keep_clust_df = sub_df[sub_df["cluster"].isin(keep_clust)]
        keep_clust_df.loc[:, "cluster"] = keep_clust_df.loc[:, "cluster"].astype("category")

        # Plot seaborn scatterplot
        plt.figure(figsize=(10, 6))
        sns.scatterplot(keep_clust_df, x=X, y=Y, hue="cluster", linewidth=0, **kwargs)
        sns.scatterplot(
            sub_df[~sub_df["cluster"].isin(keep_clust)],
            alpha=0.5,
            x=X,
            y=Y,
            linewidth=0,
            color="gray",
            **kwargs,
        )

        # Calculer les limites des axes
        max_range = max(max(sub_df[X]), max(sub_df[Y]))
        min_range = min(min(sub_df[X]), min(sub_df[Y]))
        range_val = max(abs(max_range), abs(min_range))
        buffer = 0.1 * range_val  # 10% de la plage des valeurs

        ## Définir les limites des axes avec le tampon
        # plt.xlim(-range_val - buffer, range_val + buffer)
        # plt.ylim(-range_val - buffer, range_val + buffer)

        ## Ajouter des titres et des labels
        plt.axhline(0, color="black", linestyle="--", linewidth=1)
        # Ligne horizontale à y=0
        plt.axvline(0, color="black", linestyle="--", linewidth=1)
        # Ligne verticale à x=0
        plt.title(f"{X} vs {Y} - {pdb}")
        plt.xlabel(X)
        plt.ylabel(Y)

        ## Afficher la légende avec un titre plus lisible
        if show_legend:
            plt.legend(title="Clusters", loc="upper right")
        else:
            plt.legend([], [], frameon=False)

        ## Afficher le graphique
        plt.grid(True)
        plt.show()
