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

Since the publication of Alphafold 2 software in 2021 [@Jjumper_highly_2021], the scientific community can
access a prediction accuracy of protein structures that was previously unattainable. The Alphafold 2
software has been widely adopted by the scientific community. The derivative Alphafold Multimer [@Evans2021.10.04.463034]
has been developed to predict the structure of protein complexes, once again defining a new standard for
protein-protein and protein-peptide docking. However as shown by Björn Wallner, it is sometime necessary
to generate thousands of model to get a good model [@10.1093/bioinformatics/btad573]. The analysis of the results
can then be a challenge for many researchers. Here, we present Af2-analysis, a Python package that provides
tools for the analysis of Alphafold 2 results. Af2-analysis is designed to facilitate the analysis of
protein structures predicted by Alphafold 2, providing functions for the comparison of predicted
structures with experimental structures, the visualization of predicted structures, and the
calculation of structural quality metrics. Af2-analysis is open-source and freely available on GitHub.

# Statement of need

`Af2-analysis` is a Python package that provides tools for the analysis of protein structures predicted by Alphafold 2.
The package is designed to facilitate the analysis of protein structures predicted by Alphafold 2, providing functions
for the comparison of predicted structures with experimental structures, the visualization of predicted structures, and the
calculation of structural quality metrics. The package is open-source and freely available on GitHub.


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References