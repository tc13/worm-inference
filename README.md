helminth inference model
========================

Code and data to reproduce analysis in "Diagnosis of helminths is dependent on worm fecundity and the distribution of parasites within hosts" by Crellen et al. be run using R. 
Models the relationship between worm burden and egg output as a power law (PL) function with seperate relationships for surveys in Thailand and Laos.

Code for the inference model is contained within `helminth-inference.stan`, which is run using CmdStanR from `run-helminth-model.R`

The posterior model output in .rds format can then be parsed using Rscripts `Figure1.R` and `Figure2.R` 
