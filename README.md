[![Documentation Status](https://readthedocs.org/projects/af2-analysis/badge/?version=latest)](https://af2-analysis.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/samuelmurail/af2_analysis/graph/badge.svg?token=WOJYQKKOP7)](https://codecov.io/gh/samuelmurail/af2_analysis)
[![Build Status](https://dev.azure.com/samuelmurailRPBS/af2_analysis/_apis/build/status%2Fsamuelmurail.af2_analysis?branchName=main)](https://dev.azure.com/samuelmurailRPBS/af2_analysis/_build/latest?definitionId=2&branchName=main)
[![PyPI - Version](https://img.shields.io/pypi/v/af2-analysis)](https://pypi.org/project/af2-analysis/)
[![Downloads](https://static.pepy.tech/badge/af2-analysis)](https://pepy.tech/project/af2-analysis)

# Alphafold2 Analysis

<img src="https://raw.githubusercontent.com/samuelmurail/af2_analysis/master/docs/source/logo.jpeg" alt="AF2 Analysis Logo" width="200" style="display: block; margin: auto;"/>

`af2_analysis` is a python package allowing a simplified analysis of [Alphafold][1] and [Colabfold][2] results.


## Installation

```bash
git clone https://github.com/samuelmurail/af2_analysis
cd af2_analysis
python setup.py install
```


## Usage


## Importing data

Create the `Data` object, giving the path of the directory containing the results of the alphafold2/colabfold run. 

```python
import af2_analysis
my_data = af2_analysis.Data('MY_AF2_RESULTS_DIR')
```

Extracted data are available in the `df` attribute of the `Data` object. 

```python
my_data.df
```

## Analysis

- The analysis package contains several function to add metrics like [pdockQ][3] and [pdockQ2][4]:

```python
from af2_analysis import analysis
analysis.pdockq(my_data)
analysis.pdockq(my_data)
```

## Plots

- plot msa

```python
my_data.plot_msa()
```

- plot plddt:

```python
my_data.plot_plddt([0,1])
```

- plot PAE:

```python
my_data.plot_pae(my_data.df['ranking_confidence'].idxmax())
```

- show 3D structure (`nglview` required):

```python
my_data.show_3d(my_data.df['ranking_confidence'].idxmax())
```

# References

- Jumper et al. Nature Methods (2021) doi: [10.1038/s41586-021-03819-2][1]
- Mirdita et al. Nature (2022) doi: [10.1038/s41592-022-01488-1][2]
- Bryant et al. Nat. Commun. (2022) doi: [10.1038/s41467-022-29480-5][3]
- Zhu et al. Bioinformatics (2023) doi: [10.1093/bioinformatics/btad424][4]


[1]: https://www.nature.com/articles/s41586-021-03819-2 "Jumper et al. Nature Methods (2021) doi: 10.1038/s41586-021-03819-2"
[2]: https://www.nature.com/articles/s41592-022-01488-1 "Mirdita et al. Nature (2022) doi: 10.1038/s41592-022-01488-1"
[3]: https://www.nature.com/articles/s41467-022-29480-5#citeas "Bryant et al. Nat. Commun. (2022) doi: 10.1038/s41467-022-29480-5"
[4]: https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714 "Zhu et al. Bioinformatics (2023) doi: 10.1093/bioinformatics/btad424"

