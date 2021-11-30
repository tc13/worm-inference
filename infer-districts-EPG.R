#Run FECT-Antigen worm inference model on all district
#Infers worm-egg relationship from published data
#thomas.crellen@bdi.ox.ac.uk, July 2021

#Clear workspace (if necessary)
rm(list = ls(all.names = TRUE))

#Set working directory with hard link
#setwd("/gpfs3/users/hollingsworth/mxh325/worm_inference/")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#rstan has specific installation requirements - https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(cmdstanr)
library(reshape2)
library(posterior)

#Read in data
d <- read.csv("/Users/tomc/Google Drive/Opisthorchis/Baseline/version-3/OV_cleaned_data_V3.csv") #this dataset has been updated

#Remove observations where observations are missing
d.cc <- d[!is.na(d$village_id)&!is.na(d$OV_EPG)&!is.na(d$OV_antigen_result),]

#Sort data by id and village
d.sort <- d.cc[order(d.cc$id, d.cc$village_id),]

#antigen indicator
antigen_cutoff <- ifelse(d.sort$OV_antigen_conc>=19.4,1,0) #published cut-off

#Format slides for model
slides <- rep(1, nrow(d.sort))
slide1 <- head(c(1, (cumsum(slides)+1)), -1)
slide2 <- cumsum(slides)

#Read in worm fecundity data
wf <- read.csv("/Users/tomc/Google Drive/Opisthorchis/Spencer/worm-fecund-aut-lao.csv")

#Data list for stan
in.nb <- list(N=nrow(d.sort),
              K=nrow(d.sort),
              N_clusters=length(unique(d.sort$district_id)),
              cluster=d.sort$district_id,
              slides=slides,
              slide1=slide1,
              slide2=slide2,
              eggs=d.sort$OV_EPG_int,
              urine=rep(1, nrow(d.sort)),
              antigen=antigen_cutoff,
              big_int=500,
              egg_neg=d.sort$OV_FECT_Result,
              stool_mass=rep(1,nrow(d.sort)),
              #worm fecundity data
              Nw=nrow(wf),
              E=wf$epg,
              W=wf$worms,
              N_studies=max(wf$study_id),
              study=wf$study_id)

#Total random effects in model
reffects <- in.nb$N_clusters + in.nb$N_studies

#Compile
mod1 <- cmdstan_model("/Users/tomc/Documents/worm-inference-model/worm-infer-PL.stan")

#Initial values
init_fun <- function(chain_id) list(
  #M=c(rep(1, in.nb$N_clusters), rep(147, in.nb$N_studies)),
  M_alpha=2, M_beta=0.15,
  k=0.1, y1=5, gamma=0.8,
  sp=0.91, se=0.86)

#Run HMC chains
fit <- mod1$sample(data=in.nb, refresh=10, 
                   chains=1, parallel_chains=1,
                   thin=1, seed=123, init=init_fun,
                   iter_warmup=25,
                   iter_sampling=100)

#Print summary
fit

#Extract array
post_array <- fit$draws()
fit$df <- as_draws_df(draws_array)
