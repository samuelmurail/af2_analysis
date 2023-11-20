# Basic example

Create the `Data` object, giving the path of the directory containing the results of the alphafold2/colabfold run. 

```python
import af2_analysis
my_data = af2_analysis.Data('MY_AF2_RESULTS_DIR')
```

Extracted data are available in the `df` attribute of the `Data` object. 

```python
my_data.df
```

- Compute pdockQ and pdockQ2:

 ```python
my_data.compute_pdockq()
my_data.compute_pdockq2()
 ```

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
