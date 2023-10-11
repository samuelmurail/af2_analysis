#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import pandas as pd
import json

import seaborn as sns
import matplotlib.pyplot as plt
from cmcrameri import cm

from .format import colabfold_1_5, default

class Data:
    """ Data class

    Parameters
    ----------
    dir : str
        Path to the directory containing the `log.txt` file.
    format : str
        Format of the data.
    df : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    """
    
    def __init__(self, directory=None):
        """
        """

        if directory is not None:
            self.read_directory(directory)
    
    def read_directory(self, directory, keep_recycles=False):
        """ Read a directory.

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

        if os.path.isfile(os.path.join(directory, 'log.txt')):
            self.format = 'colabfold_1.5'
            self.df = colabfold_1_5.read_log(directory, keep_recycles)
            self.add_pdb()
            self.add_json()
        else:
            self.format = 'default'
            self.df = default.read_dir(directory)
            self.add_json()
            self.extract_json()

    def add_json(self):
        """ Add json files to the dataframe.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == 'colabfold_1.5':
            colabfold_1_5.add_json(self.df, self.dir)
        if self.format == 'default':
            default.add_json(self.df, self.dir)
    
    def extract_json(self):
        """ Extract json files to the dataframe.
        
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
                
                with open(row.json, 'r') as f:
                    data = json.load(f)
                
                json_list.append(data)
        
        new_column = {}
        for keys in json_list[0].keys():
            new_column[keys] = []
        for data in json_list:
            for keys in data.keys():
                new_column[keys].append(data[keys])
        
        for keys in new_column.keys():
            self.df[keys] = np.nan
            new_col = pd.Series(new_column[keys], index=index_list)
            self.df[keys].iloc[index_list] = new_col

    def add_pdb(self):
        """ Add pdb files to the dataframe.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.format == 'colabfold_1.5':
            colabfold_1_5.add_pdb(self.df, self.dir)

    def add_fasta(self, csv):
        """ Add fasta sequence to the dataframe.
        
        Parameters
        ----------
        csv : str
            Path to the csv file containing the fasta sequence.

        Returns
        -------
        None
        """
            
        if self.format == 'colabfold_1.5':
            colabfold_1_5.add_fasta(self.df, csv)
    
    def keep_last_recycle(self):
        """ Keep only the last recycle for each query.

        """

        idx = self.df.groupby(['query', 'seed', 'model', 'weight'])['recycle'].transform(max) == self.df['recycle']
        self.df = self.df[idx]


    def plot_maxscore_as_col(self, score, col, hue='query'):
        
        col_list = self.df[col].unique()
        query_list = self.df[hue].unique()
        #print(col_list)
        #print(query_list)
        
        out_list = []
        
        for query in query_list:
            #print(query)
            query_pd = self.df[ self.df[hue] == query]
            
            for column in col_list:
                #print(column)
                #~print()
                
                col_pd = query_pd [ query_pd[col] <= column ]
                #print(col_pd[score])
                #print(column, len(col_pd))
                #print(col, col_pd.columns)
                
                if len(col_pd) > 0:
                
                    out_list.append({
                        hue: query,
                        score: col_pd[score].max(),
                        col: column})
                    #print(column, len(col_pd), col_pd[score].max())
        
        max_pd = pd.DataFrame(out_list)
    
        ax = sns.lineplot(max_pd, x=col, y=score, hue=hue)

        return(ax)
    
    def plot_pae(self, index, cmap=cm.vik):

        row = self.df.iloc[index]

        if row['json'] is None:
            return(None, None)

        with open(row['json']) as f:
            local_json = json.load(f)
        
        pae_array = np.array(local_json['pae'])

        fig = plt.figure()
        ax = plt.imshow(pae_array, cmap=cmap, vmin=0., vmax=30.)
        plt.colorbar(ax)

        return(fig, ax)