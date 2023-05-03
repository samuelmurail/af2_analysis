#!/usr/bin/env python3
# coding: utf-8

import os
import re
import logging
from tqdm.auto import tqdm
import pandas as pd

logger = logging.getLogger()
logger.setLevel(logging.INFO)

weigths = [
    'alphafold2',
    'alphafold2_ptm',
    'alphafold2_multimer_v1',
    'alphafold2_multimer_v2',
    'alphafold2_multimer_v3'
]

def read_log(directory):
    """ Read colabfold `log.txt` file. 
    Extract information about the model and the query.

    Return a pandas dataframe.

    Parameters
    ----------
    directory : str
        Path to the directory containing the `log.txt` file.
    
    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.

    """
    
    logger.info(f"Reading {os.path.join(directory, 'log.txt')}")

    log_file = os.path.join(directory, "log.txt")
    
    with open(log_file) as f:
        log_lines = f.readlines()
    
    log_dict_list = []
    
    for line in log_lines:
        if line.find("Query") != -1:
            query = line.split()[4]
            #print(query)
        if line.find("recycle=") != -1:

            token = line.split()
            
            name_token = token[2].split("_")
            weight = "_".join(name_token[:-4])
            assert weight in weigths
            model = int(name_token[-3])
            seed = int(name_token[-1])
            recycle = int(token[3].split("=")[1])
            pLDDT= float(token[4].split("=")[1])
            pTM= float(token[5].split("=")[1])
            ipTM= float(token[6].split("=")[1])
            
            log_dict_list.append({
                'query': query,
                'seed': seed,
                'model': model,
                'weight': weight,
                'recycle': recycle,
                'pLDDT': pLDDT,
                'pTM': pTM,
                'ipTM': ipTM,
            })
    
    log_pd = pd.DataFrame(log_dict_list)
    log_pd['ranking'] = 0.8 * log_pd['ipTM'] + 0.2 * log_pd['pTM'] 
    return log_pd


def add_pdb(log_pd, directory):
    """ Find pdb files in the directory.

    Parameters
    ----------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    directory : str
        Path to the directory containing the pdb files.
    
    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    #logger.info(f"Extracting pdb files location")
    raw_list = os.listdir(directory)
    file_list = []
    for file in raw_list:
        if file.endswith(".pdb"):
            file_list.append(file)
    
    pdb_list = []

    for _, row in tqdm(log_pd.iterrows(), total=log_pd.shape[0]):
                
        reg = fr"{row['query']}_.*_{row['weight']}_model_{row['model']}_seed_{row['seed']:03d}\.r{row['recycle']}\.pdb"
        r = re.compile(reg)
        res = list(filter(r.match, file_list))

        if len(res) == 1:
            pdb_list.append(os.path.join(directory, res[0]))
            file_list.remove(res[0])
        elif len(res) > 1:
            logger.warning(f"Multiple pdb file for {reg}: {res}")
            pdb_list.append(None)
        else:
            logger.info(f"Not founded : {reg}")
            pdb_list.append(None)
    
    log_pd.loc[:, 'pdb'] = pdb_list

def add_json(log_pd, directory):
    """ Find json files in the directory.

    Parameters
    ----------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    directory : str
        Path to the directory containing the json files.
    
    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    logger.info(f"Extracting json files location")

    raw_list = os.listdir(directory)
    file_list = []
    for file in raw_list:
        if file.endswith(".json"):
            file_list.append(file)
    
    json_list = []
    
    # Get the last recycle:
    last_recycle = log_pd.groupby(['query', 'seed', 'model', 'weight'])['recycle'].transform(max) == log_pd['recycle']
    
    for i, last in tqdm(enumerate(last_recycle), total=len(last_recycle)):
        row = log_pd.iloc[i]
        file_state = False
        
        if not last:
            json_list.append(None)
            continue
        
        reg = fr"{row['query']}_scores_.*_{row['weight']}_model_{row['model']}_seed_{row['seed']:03d}\.json"
        r = re.compile(reg)
        
        res = list(filter(r.match, file_list))

        if len(res) == 1:
            json_list.append(os.path.join(directory, res[0]))
            #file_list.remove(res[0])
        elif len(res) == 0:
            logger.warning(f"Not founded : {reg}")
            json_list.append(None)
        else:
            logger.warning(f"Multiple json file for {reg}: {res}")
            json_list.append(res[0])
            #file_list.remove(res[0])
    
    log_pd.loc[:, 'json'] = json_list

def add_fasta(log_pd, csv):
    """ Add the raw sequence from csv file in the dataframe.

    Parameters
    ----------
    log_pd : pandas.DataFrame
        Dataframe containing the information extracted from the `log.txt` file.
    csv : str
        Path to the csv file containing the raw sequence.

    Returns
    -------
    None
        The `log_pd` dataframe is modified in place.
    """

    logger.info(f"Extracting sequence from csv file")

    seq_pd = pd.read_csv(csv)
    seq_dict = {}

    for _, seq_row in seq_pd.iterrows():
        seq_dict[seq_row['id']] = seq_row['sequence']

    seq_list = []
    for _, row in tqdm(log_pd.iterrows(), total=log_pd.shape[0]):
        seq_list.append(seq_dict[row['query']])

    log_pd.loc[:, 'sequence'] = seq_list

