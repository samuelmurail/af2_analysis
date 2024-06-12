import numpy as np

from tqdm.auto import tqdm
from . import data

def extract_pae_pep(my_data, fun = np.mean):
    pep_rec_pae_list = []
    rec_pep_pae_list = []
    
    for i, (query, json) in tqdm(
        enumerate(zip(my_data.df['query'], my_data.df['json'])), total=len(my_data.df)):
        
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        pae = data.get_pae(json)
        
        #print(f"{cum_sum_chain}  {cum_sum_chain[-2]}:{cum_sum_chain[-1]}")
        rec_pep_pae = fun(pae[0:cum_sum_chain[-2], cum_sum_chain[-2]:cum_sum_chain[-1]])
        pep_rec_pae = fun(pae[cum_sum_chain[-2]:cum_sum_chain[-1], 0:cum_sum_chain[-2]])
        
        pep_rec_pae_list.append(pep_rec_pae)
        rec_pep_pae_list.append(rec_pep_pae)
    
    my_data.df['PAE_pep_rec'] = pep_rec_pae_list
    my_data.df['PAE_rec_pep'] = rec_pep_pae_list


def extract_plddt_pep(my_data, fun = np.mean):
    pep_plddt_list = []
    
    for i, (query, pdb) in tqdm(enumerate(zip(my_data.df['query'], my_data.df['pdb'])), total=len(my_data.df)):
        
        chain_length = my_data.chain_length[query]
        cum_sum_chain = np.cumsum([0] + chain_length)

        plddt = my_data.get_plddt(i)
        pep_plddt_list.append(fun(plddt[cum_sum_chain[-2]:cum_sum_chain[-1]]))
    
    my_data.df['plddt_pep'] = pep_plddt_list

def compute_LIS_pep(my_data):

    my_data.compute_LIS_matrix()

    pep_LIS_list = []
    pep_LIS2_list = []

    fun = np.max
    
    for query, LIS in zip(my_data.df['query'], my_data.df['LIS']):
        chain_num = len(my_data.chains[query])

        pep_LIS_list.append(fun(LIS[0:chain_num-1,chain_num - 1]))
        pep_LIS2_list.append(fun(LIS[chain_num - 1, 0:chain_num-1]))
    
    my_data.df['LIS_rec_pep'] = pep_LIS2_list
    my_data.df['LIS_pep_rec'] = pep_LIS_list