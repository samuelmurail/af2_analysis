#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pandas as pd
import json
import pdb_numpy
import seaborn as sns
import matplotlib.pyplot as plt
from cmcrameri import cm
from tqdm.auto import tqdm
from scipy.spatial import distance_matrix

import ipywidgets as widgets

from .format import colabfold_1_5, default
from . import sequence, plot


class Data:
    """Data class

    Parameters
    ----------
    dir : str
        Path to the directory containing the `log.txt` file.
    format : str
        Format of the data.
    df : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    chains : dict
        Dictionary containing the chains of each query.
    chain_length : dict
        Dictionary containing the length of each chain of each query.
    
    Methods
    -------
    read_directory(directory, keep_recycles=False)
        Read a directory.
    export_csv(path)
        Export the dataframe to a csv file.
    import_csv(path)
        Import a csv file to the dataframe.
    add_json()
        Add json files to the dataframe.
    extract_json()
        Extract json files to the dataframe.
    add_pdb()
        Add pdb files to the dataframe.
    add_fasta(csv)
        Add fasta sequence to the dataframe.
    keep_last_recycle()
        Keep only the last recycle for each query.
    plot_maxscore_as_col(score, col, hue='query')
        Plot the maxscore as a function of a column.
    plot_pae(index, cmap=cm.vik)
        Plot the PAE matrix.
    plot_plddt(index_list)
        Plot the pLDDT.
    show_3d(index)
        Show the 3D structure.
    compute_pdockq()
        Compute pdockq from the pdb file.
    compute_mpdockq()
        Compute mpdockq from the pdb file.
    compute_pdockq2()
        Compute pdockq2 from the pdb file.
    plot_msa(filter_qid=0.15, filter_cov=0.4)
        Plot the msa from the a3m file.
    show_plot_info()
        Show the plot info.
    extract_inter_chain_pae(fun=np.mean)
        Read the PAE matrix and extract the average inter chain PAE.

    """

    def __init__(self, directory=None, csv=None):
        """ """

        if directory is not None:
            self.read_directory(directory)
        elif csv is not None:
            self.import_csv(csv)

    def read_directory(self, directory, keep_recycles=False):
        """Read a directory.

        If the directory contains a `log.txt` file, the format is set to `colabfold_1.5`.

        Parameters
        ----------
        directory : str
            Path to the directory containing the `log.txt` file.

        Returns
        -------
        None
        """
        self.dir = directory

        if os.path.isfile(os.path.join(directory, "log.txt")):
            self.format = "colabfold_1.5"
            self.df = colabfold_1_5.read_log(directory, keep_recycles)
            self.add_pdb()
            self.add_json()
        else:
            self.format = "default"
            self.df = default.read_dir(directory)
            self.add_json()
            # self.extract_json()

        self.chains = {}
        self.chain_length = {}
        for querie in self.df["query"].unique():
            # print(querie, self.df[self.df['query'] == querie])
            first_model = pdb_numpy.Coor(
                self.df[self.df["query"] == querie].iloc[0]["pdb"]
            )
            self.chains[querie] = list(np.unique(first_model.models[0].chain))
            self.chain_length[querie] = [
                len(
                    np.unique(
                        first_model.models[0].uniq_resid[
                            first_model.models[0].chain == chain
                        ]
                    )
                )
                for chain in self.chains[querie]
            ]

    def export_csv(self, path):
        """Export the dataframe to a csv file.

        Parameters
        ----------
        path : str
            Path to the csv file.

        Returns
        -------
        None
        """

        self.df.to_csv(path, index=False)

    def import_csv(self, path):
        """Import a csv file to the dataframe.

        Parameters
        ----------
        path : str
            Path to the csv file.

        Returns
        -------
        None
        """

        self.df = pd.read_csv(path)
        self.dir = os.path.dirname(self.df['pdb'][0])

        self.chains = {}
        self.chain_length = {}
        for querie in self.df["query"].unique():
            # print(querie, self.df[self.df['query'] == querie])
            first_model = pdb_numpy.Coor(
                self.df[self.df["query"] == querie].iloc[0]["pdb"]
            )
            self.chains[querie] = list(np.unique(first_model.models[0].chain))
            self.chain_length[querie] = [
                len(
                    np.unique(
                        first_model.models[0].uniq_resid[
                            first_model.models[0].chain == chain
                        ]
                    )
                )
                for chain in self.chains[querie]
            ]

    def add_json(self):
        """Add json files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_json(self.df, self.dir)
        if self.format == "default":
            default.add_json(self.df, self.dir)

    def extract_json(self):
        """Extract json files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        index_list = []
        json_list = []

        for i, row in self.df.iterrows():
            if row.json is not None:
                index_list.append(i)

                with open(row.json, "r") as f:
                    data = json.load(f)

                json_list.append(data)

        new_column = {}
        for keys in json_list[0].keys():
            new_column[keys] = []
        for data in json_list:
            for keys in data.keys():
                new_column[keys].append(data[keys])

        for keys in new_column.keys():
            self.df[keys] = None
            new_col = pd.Series(new_column[keys], index=index_list)
            self.df[keys].iloc[index_list] = new_col

    def add_pdb(self):
        """Add pdb files to the dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_pdb(self.df, self.dir)

    def add_fasta(self, csv):
        """Add fasta sequence to the dataframe.

        Parameters
        ----------
        csv : str
            Path to the csv file containing the fasta sequence.

        Returns
        -------
        None
        """

        if self.format == "colabfold_1.5":
            colabfold_1_5.add_fasta(self.df, csv)

    def keep_last_recycle(self):
        """Keep only the last recycle for each query."""

        idx = (
            self.df.groupby(["query", "seed", "model", "weight"])["recycle"].transform(
                "max"
            )
            == self.df["recycle"]
        )
        self.df = self.df[idx]

    def plot_maxscore_as_col(self, score, col, hue="query"):
        col_list = self.df[col].unique()
        query_list = self.df[hue].unique()
        # print(col_list)
        # print(query_list)

        out_list = []

        for query in query_list:
            # print(query)
            query_pd = self.df[self.df[hue] == query]

            for column in col_list:
                # print(column)
                # ~print()

                col_pd = query_pd[query_pd[col] <= column]
                # print(col_pd[score])
                # print(column, len(col_pd))
                # print(col, col_pd.columns)

                if len(col_pd) > 0:
                    out_list.append(
                        {hue: query, score: col_pd[score].max(), col: column}
                    )
                    # print(column, len(col_pd), col_pd[score].max())

        max_pd = pd.DataFrame(out_list)

        fig, ax = plt.subplots()
        sns.lineplot(max_pd, x=col, y=score, hue=hue)

        return (fig, ax)

    def plot_pae(self, index, cmap=cm.vik):
        row = self.df.iloc[index]

        if row["json"] is None:
            return (None, None)

        with open(row["json"]) as f:
            local_json = json.load(f)

        pae_array = np.array(local_json["pae"])

        fig, ax = plt.subplots()
        res_max = sum(self.chain_length[row["query"]])

        img = ax.imshow(
            pae_array,
            cmap=cmap,
            vmin=0.0,
            vmax=30.0,
        )

        plt.hlines(
            np.cumsum(self.chain_length[row["query"]][:-1]) - 0.5,
            xmin=-0.5,
            xmax=res_max,
            colors="black",
        )

        plt.vlines(
            np.cumsum(self.chain_length[row["query"]][:-1])- 0.5,
            ymin=-0.5,
            ymax=res_max,
            colors="black",
        )

        plt.xlim(-0.5, res_max-0.5)
        plt.ylim(res_max-0.5, -0.5)
        ax.set_yticklabels(self.chains[row["query"]])
        chain_pos = []
        len_sum = 0
        for longueur in self.chain_length[row["query"]]:
            chain_pos.append(len_sum + longueur / 2)
            len_sum += longueur

        ax.set_yticks(chain_pos)
        cbar = plt.colorbar(img)
        cbar.set_label("Predicted Aligned Error (Å)", rotation=270)
        cbar.ax.get_yaxis().labelpad = 15

        return (fig, ax)

    def plot_plddt(self, index_list):
        fig, ax = plt.subplots()

        for index in index_list:
            row = self.df.iloc[index]

            if row["json"] is None:
                return (None, None)

            with open(row["json"]) as f:
                local_json = json.load(f)

            plddt_array = np.array(local_json["plddt"])

            plt.plot(plddt_array)

        plt.vlines(
            np.cumsum(self.chain_length[row["query"]][:-1]),
            ymin=0,
            ymax=100.0,
            colors="black",
        )
        plt.ylim(0, 100)
        plt.xlim(0, sum(self.chain_length[row["query"]]))
        plt.xlabel("Residue")
        plt.ylabel("predicted LDDT")

        return (fig, ax)

    def show_3d(self, index):
        row = self.df.iloc[index]

        if row["pdb"] is None:
            return (None, None)

        import nglview as nv

        # Bug with show_file
        # view = nv.show_file(row['pdb'])
        view = nv.show_structure_file(row["pdb"])
        # view.add_component(ref_coor[0])
        # view.clear_representations(1)
        # view[1].add_cartoon(selection="protein", color='blue')
        # view[1].add_licorice(selection=":A", color='blue')
        # view[0].add_licorice(selection=":A")
        return view

    def compute_pdockq(self):
        """
        Compute pdockq from the pdb file.

        """

        from pdb_numpy.analysis import compute_pdockQ

        pdockq_list = []

        for pdb in tqdm(self.df["pdb"], total=len(self.df["pdb"])):
            if pdb:
                model = pdb_numpy.Coor(pdb)
                pdockq_list += compute_pdockQ(model)
            else:
                pdockq_list.append(None)

        self.df["pdockq"] = pdockq_list


    def compute_mpdockq(self):
        """
        Compute mpdockq from the pdb file.

        """

        from pdb_numpy.analysis import compute_pdockQ

        pdockq_list = []

        for pdb in tqdm(self.df["pdb"], total=len(self.df["pdb"])):
            if (pdb is not None and pdb is not np.nan):
                model = pdb_numpy.Coor(pdb)
                pdockq_list += compute_pdockQ(
                    model,
                    cutoff=8.0,
                    L=0.728,
                    x0=309.375,
                    k=0.098,
                    b=0.262)
            else:
                pdockq_list.append(None)

        self.df["mpdockq"] = pdockq_list

    def compute_pdockq2(self):
        """
        Compute pdockq2 from the pdb file.

        $$ pDockQ_2 = \frac{L}{1 + exp [-k*(X_i-X_0)]} + b$$

        with

        $$ X_i = \langle \frac{1}{1+(\frac{PAE_{int}}{d_0})^2} \rangle - \langle pLDDT \rangle_{int}$$

        Ref:
        https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
        """

        from pdb_numpy.analysis import compute_pdockQ2

        pdockq_list = []

        max_chain_num = 0
        for query in self.chains:
            chain_num = len(self.chains[query])
            if chain_num > max_chain_num:
                max_chain_num = chain_num

        for i in range(max_chain_num):
            pdockq_list.append([])

        for pdb, json_path in tqdm(
            zip(self.df["pdb"], self.df["json"]), total=len(self.df["pdb"])
        ):
            if (pdb is not None and pdb is not np.nan and json_path is not None and json_path is not np.nan):
                model = pdb_numpy.Coor(pdb)
                with open(json_path) as f:
                    local_json = json.load(f)
                pae_array = np.array(local_json["pae"])

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
            self.df[f"pdockq2_{chr(65+i)}"] = pdockq_list[i]


    def compute_piTM(self):
        r"""Compute the piTM score as define in [2]_.

        .. math::
            piTM = \max_{i \in \mathcal{I}} \frac{1}{I} \sum_{j \in \mathcal{I}}  \frac{1}{1 + [\langle e_{ij} \rangle / d_0 (I)]^2}

        with:

        .. math::
            d_0(I) = \begin{cases} 1.25 \sqrt[3]{I -15} -1.8\text{,} & \text{if } I \geq 22 \\ 0.02 I \text{,} & \text{if } I < 22  \end{cases}
        

        Implementation was inspired from `predicted_tm_score_v1()` in https://github.com/FreshAirTonight/af2complex/blob/main/src/alphafold/common/confidence.py

        References
        ----------
        .. [2] Mu Gao, Davi Nakajima An, Jerry M. Parks & Jeffrey Skolnick. 
            AF2Complex predicts direct physical interactions in multimeric proteins with deep learning
            *Nature Communications*. volume 13, Article number: 1744 (2022).
            https://www.nature.com/articles/s41467-022-29394-2
        """

        from pdb_numpy.analysis import compute_piTM

        piTM_chain_list = []
        piTM_list = []

        max_chain_num = 0
        for query in self.chains:
            chain_num = len(self.chains[query])
            if chain_num > max_chain_num:
                max_chain_num = chain_num

        for i in range(max_chain_num):
            piTM_chain_list.append([])

        for pdb, json_path in tqdm(
            zip(self.df["pdb"], self.df["json"]), total=len(self.df["pdb"])
        ):
            # print(pdb, json_path)
            if (pdb is not None and pdb is not np.nan and json_path is not None and json_path is not np.nan):
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
        self.df["piTM"] = piTM_list
        for i in range(max_chain_num):
            self.df[f"piTM_{chr(65+i)}"] = piTM_chain_list[i]


    def plot_msa(self, filter_qid=0.15, filter_cov=0.4):
        """
        Plot the msa from the a3m file.

        Parameters
        ----------
        filter_qid : float
            Minimal sequence identity to keep a sequence.
        filter_cov : float
            Minimal coverage to keep a sequence.

        Returns
        -------
        None

        ..Warning only tested with colabfold 1.5
        """

        raw_list = os.listdir(self.dir)
        file_list = []
        for file in raw_list:
            if file.endswith(".a3m"):
                file_list.append(file)

        for a3m_file in file_list:
            print(f"Reading MSA file:{a3m_file}")
            querie = a3m_file.split("/")[-1].split(".")[0]

            a3m_lines = open(os.path.join(self.dir, a3m_file), "r").readlines()[1:]
            seqs, mtx, nams = sequence.parse_a3m(
                a3m_lines=a3m_lines,
                filter_qid=filter_qid,
                filter_cov=filter_cov)
            print(f"- Keeping {len(seqs):6} sequences for plotting.")
            feature_dict = {}
            feature_dict["msa"] = sequence.convert_aa_msa(seqs)
            feature_dict["num_alignments"] = len(seqs)

            if len(seqs) == sum(self.chain_length[querie]):
                feature_dict["asym_id"] = []
                for i, chain_len in enumerate(self.chain_length[querie]):
                    feature_dict["asym_id"] += [i + 1.] * chain_len
                feature_dict["asym_id"] = np.array(feature_dict["asym_id"])

            fig = plot.plot_msa_v2(feature_dict)
            plt.show()

    def show_plot_info(self):
        """
        Need to solve the issue with:

        ```
        %matplotlib ipympl
        ```

        plots don´t update when changing the model number.

        """

        model_widget = widgets.IntSlider(
            value=1,
            min=1,
            max=len(self.df),
            step=1,
            description='model:',
            disabled=False,
        )
        display(model_widget)

        def show_model(rank_num):
            print(rank_num)
            plddt_fig, plddt_ax = self.plot_plddt([rank_num - 1])
            pae_fig, pae_ax = self.plot_pae(rank_num - 1)
            view = self.show_3d(rank_num - 1)
            #view.zoomTo()
            return view
            
        output = widgets.Output()
        display(output)
        with output:
            show_model(model_widget.value)
            #logger.info(results['metric'][0][rank_num - 1]['print_line'])

        def on_value_change(change):
            output.clear_output()
            with output:
                show_model(model_widget.value)

        model_widget.observe(on_value_change, names='value')


    def extract_inter_chain_pae(self, fun=np.mean):
        """ Read the PAE matrix and extract the average inter chain PAE.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        pae_list = []
        
        for query, json_path in tqdm(
                        zip(self.df["query"], self.df["json"]), total=len(self.df["json"])
            ):
            if (json_path is not None and json_path is not np.nan):
                with open(json_path) as f:
                    local_json = json.load(f)
                pae_array = np.array(local_json["pae"])
                
                chain_lens = self.chain_length[query]
                chain_len_sums = np.cumsum([0] + chain_lens)
                pae_chain_array = np.empty((len(chain_lens), len(chain_lens)))
                chain_ids = self.chains[query]
                
                pae_dict = {}
                
                for i in range(len(chain_lens)):
                    for j in range(len(chain_lens)):
                        pae_val = fun(pae_array[chain_len_sums[i]:chain_len_sums[i+1],chain_len_sums[j]:chain_len_sums[j+1]])
                        pae_dict[f"PAE_{chain_ids[i]}_{chain_ids[j]}"] = pae_val
                
                pae_list.append(pae_dict)
            else:
                pae_list.append({})

        pae_df = pd.DataFrame(pae_list)
        
        for col in pae_df.columns:
            self.df[col] = pae_df[col]


def concat_data(data_list):
    """ Concatenate data from a list of Data objects.

    Parameters
    ----------
    data_list : list
        List of Data objects.

    Returns
    -------
    Data
        Concatenated Data object.   
    """

    concat = Data(directory=None, csv=None)

    concat.df = pd.concat([data.df for data in data_list], ignore_index=True)
    concat.chains = data_list[0].chains
    concat.chain_length = data_list[0].chain_length
    for i in range(1, len(data_list)):
        concat.chains.update(data_list[i].chains)
        concat.chain_length.update(data_list[i].chain_length)
    
    return concat