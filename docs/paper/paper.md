---
title: 'Af2-analysis: a Python package for alphafold 2 analysis'
tags:
  - Python
  - Alphafold
  - Protein Structure
  - Structural bioinformatics
authors:
  - name: Alaa Reguei
    affiliation: "1"
  - name: Samuel Murail
    orcid: "0000-0002-3842-5333"
    affiliation: "1"

affiliations:
 - name: "Université Paris Cité, Inserm U1133, CNRS UMR 8251, Paris, France"
   index: 1

date: 11 July 2024
bibliography: paper.bib
---

# Summary

The publication of AlphaFold 2[@jumper_highly_2021] has significantly advanced the field of protein structure prediction. Predicting protein structures has long been a central challenge in structural bioinformatics, with the ultimate aim of elucidating the relationship between protein structure and function. Accurate protein structure prediction is crucial for various applications, including drug discovery, protein engineering, and the investigation of protein-protein interactions. AlphaFold 2, a deep learning-based approach, has demonstrated unprecedented accuracy in protein structure prediction, surpassing other contemporary methods. Here, we present `af2-analysis`, a Python package that provides tools for the analysis of Alphafold 2 results. `af2-analysis` is designed to facilitate the analysis of protein structures predicted by Alphafold 2, providing functions for the comparison of predicted structures with experimental structures, the visualization of predicted structures, and the calculation of structural quality metrics.

# Statement of need

Since the publication of Alphafold 2[@jumper_highly_2021] in 2021, the scientific community can access a prediction accuracy of protein structures that was previously unattainable. The derivatives Colabfold[@mirdita2022colabfold] and Alphafold Multimer[@Evans2021.10.04.463034] have been developed to predict the structure of protein complexes, once again defining a new standard for protein-protein and protein-peptide docking. However as shown by Björn Wallner, it is sometime necessary to generate thousands of model to get few good models[@10.1093/bioinformatics/btad573]. The analysis of the results can then be a fastidious, moreover if alphafold quality metrics are good, supplementary metrics have been developed to evaluate the quality of the models as pdockq[@bryant2022improved], pdockq2[@10.1093/bioinformatics/btad424], and LIS score[@Kim2024.02.19.580970]. All this metrics have to be calculated from different scripts. Another point to access is the diversity of the models, and as shown in AFsample[@10.1093/bioinformatics/btad573] it is sometime necessary to compute up to ten thousand models, and then to cluster them to select the best models. `Af2-analysis` is designed to facilitate the analysis of set of model structures and metrics. The library is build upon `pandas` library and will import an alphafold or colabfold prediction directory as a `pandas` dataframe. The library will provide functions to add additional metric to the dataframe, to compare the models with experimental structures, to visualize the models, to cluster the models, and to select the best models.

# References