helminth inference model
========================

Code and data to reproduce analysis in "[Diagnosis of helminths is dependent on worm fecundity and the distribution of parasites within hosts](https://royalsocietypublishing.org/doi/10.1098/rspb.2022.2204)" by Crellen T _et_ _al_. 2023 _Proc_ _Roy_ _Soc_ _B_ be run using R and Stan. 

This models the relationship between total helminth worm burden within a host and egg output as a power law (PL) function. The code can model seperate relationships for surveys in different populations (Thailand and Laos in this analysis).

Code for the inference model is contained within `helminth-inference.stan`, which is run using CmdStanR from `run-helminth-model.R`

The posterior model output in .rds format can then be parsed using Rscripts `Figure1.R` and `Figure2.R` to replicate figures in the paper.

The file containing the necessary dataset `datasets.csv` must be downloaded [from Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.q83bk3jn6). 
