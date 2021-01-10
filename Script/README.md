# Script

This subdirectory repository contains the scripts for preparing and modelling the datasets.

- [DATRAS_filter.R](DATRAS_filter.R) - contains the R code to go from the raw ICES cpue length per hour to filtered hauls.
- [PCMCI_prep.R](PCMCI_prep.R) - contains the R code to go from the filtered hauls to the aggregated time series used in the PCMCI analysis.
- [PCMCI.ipynb](PCMCI.ipynb) - contains a Jupyter notebook for running PCMCI analysis on the aggregated time series prepared in PCMCI_prep.R
- [Patterns_in_causal_association_network.ipynb](Patterns_in_causal_association_network.ipynb) - contains a Jupyter notebook to bring together data on life history traits, PCMCI results and differences in juvenile and adult biomass, that is used in Random Forest analysis.
- [Random_Forest_Shapley.ipynb](Random_Forest_Shapley.ipynb)- - Performs the Random Forest analysis and extracts feature importance based on Shapley values.

 <br>
 <br>
 
![](images/cover.png)
