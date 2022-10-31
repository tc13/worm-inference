#!/usr/bin/env

#Run autopsy / expulsion worm inference model
#Power Law function (PL), negative binomial hurdle likelihood, autopsy epg divided by factor
#Stoll factor estimated
#Two curves (TH and LAO)

#Libraries
require(cmdstanr)
require(posterior)

#Model code
mod <- cmdstan_model("helminth-inference.stan")

#Read in cleaned data
d <- read.csv("datasets.csv")

#Add variables
studies <- unique(d$survey)
d$study_id <- match(d$survey, studies)
countries <- unique(d$country)
d$country_id <- match(d$country, countries)
d$study_type_id <- ifelse(d$study_type=="autopsy",0,1)

survey.df <- data.frame(
  name=studies,
  id=match(studies, studies),
  mean_worms = sapply(studies, FUN=function(x) mean(d$worms[d$survey==x])),
  study_pop = sapply(studies, FUN=function(x) d$country_id[which(d$survey==x)][1]),
  ecc_id = ifelse(studies=="TH2",1,0)
  )

#Format list for stan input
input <- list(
  N = nrow(d),
  epg = d$epg,
  worms = d$worms,
  N_studies = max(d$study_id),
  study_id = d$study_id,
  study_type = d$study_type_id,
  mean_worms_obs = survey.df$mean_worms,
  N_pops = length(countries),
  study_pop = survey.df$study_pop,
  N_expul = sum(d$study_type_id),
  N_autopsy = nrow(d) - sum(d$study_type_id),
  expul_indx = which(d$study_type=="expulsion"),
  autopsy_indx = which(d$study_type=="autopsy"),
  ecc_id = survey.df$ecc_id,
  model_function = 1, #Power Law function
  error_dist = 3, #Negative Binomial Hurdle error distribution
  delta_worm = 2500
  )

#Run model
chains <- mod$sample(data=input, 
                     chains=4, 
                     parallel_chains=4,
                     seed=128,
                     iter_warmup=2750,
		     #iter_warmup=725,
		     adapt_delta=0.82,
                     iter_sampling=250,
                     sig_figs=5,
		     output_dir="posteriors",
		     output_basename="wm-pl-nbh-2p")

#Vector of parameter names
vars <- c("M[1]","M[2]","M[3]", "M[4]", "M[5]" ,"k[1]","k[2]", "k[3]", "k[4]","k[5]" ,"pr_recovery", "y1[1]", "gamma[1]", "y1[2]", "gamma[2]","k_mean", "k_sd", "h[1]", "h[2]" ,"b[1]", "ec_factor")

#Print MCMC diagnostics
chains$print(vars, max_rows=length(vars)+1)

#Read in chains
fit <-  read_cmdstan_csv(files=c("posteriors/wm-pl-nbh-2p-1.csv",
				 "posteriors/wm-pl-nbh-2p-2.csv",
				 "posteriors/wm-pl-nbh-2p-3.csv",
				 "posteriors/wm-pl-nbh-2p-4.csv"),
				 variables = c(vars, "aut_state", "pState", "log_lik"),
				 sampler_diagnostics = "")

#Indexing
indx.pars <- length(vars)
indx.aut <- input$N_autopsy*4
indx.expul <- input$N_expul*input$delta_worm
start.ll <- indx.pars+indx.aut+indx.expul
end.ll <- dim(fit$post_warmup_draws)[3]

#Extract main parameters as posterior data.frame
post <- as_draws_df(fit$post_warmup_draws[,,vars])

#Extract inferred states - i) true worm burden in expulsion studies, and ii) multiplication factor for EPG from autopsy study
aut_matrix  <- as_draws_matrix(fit$post_warmup_draws[,,c((indx.pars+1):(indx.pars+indx.aut))])
expul_matrix<- as_draws_matrix(fit$post_warmup_draws[,,c((indx.pars+indx.aut+1):(indx.pars+indx.aut+indx.expul))])

pState <- apply(expul_matrix, 2, mean)
aState <- apply(aut_matrix, 2, mean)

e.vec <- a.vec <- c()

for(i in 1:input$N_expul){
	indx.temp <- which(grepl(paste("[",i,",",sep=""), names(pState), fixed=TRUE))
	p.max <- which(pState[indx.temp]==max(pState[indx.temp]))
	e.vec[i] <- p.max
}

for(i in 1:input$N_autopsy){
	indx.temp <- which(grepl(paste("[",i,",",sep=""), names(aState), fixed=TRUE))
	p.max <- which(aState[indx.temp]==max(aState[indx.temp]))
	a.vec[i] <- ifelse(length(p.max)>1, 2, p.max)
}

#Extract log liklihood as array
log_lik <- as_draws_array(fit$post_warmup_draws[,,c((start.ll+1):end.ll)])

#Collect variables into a list
posterior.list <- list(
	M1 = unlist(post["M[1]"], use.names=F),
	M2 = unlist(post["M[2]"], use.names=F),
	M3 = unlist(post["M[3]"], use.names=F),
	M4 = unlist(post["M[4]"], use.names=F),
	M5 = unlist(post["M[5]"], use.names=F),
	k1 = unlist(post["k[1]"], use.names=F),
	k2 = unlist(post["k[2]"], use.names=F),
	k3 = unlist(post["k[3]"], use.names=F),
	k4 = unlist(post["k[4]"], use.names=F),
	k5 = unlist(post["k[5]"], use.names=F),
	k_mu = unlist(post["k_mean"], use.names=F),
	k_sd = unlist(post["k_sd"], use.names=F),
	pr_recovery = unlist(post["pr_recovery"], use.names=F),
	y1_TH = unlist(post["y1[1]"], use.names=F),
	gamma_TH = unlist(post["gamma[1]"], use.names=F),
	h_TH = unlist(post["h[1]"], use.names=F),
	y1_LAO = unlist(post["y1[2]"], use.names=F),
        gamma_LAO = unlist(post["gamma[2]"], use.names=F),
        h_LAO = unlist(post["h[2]"], use.names=F),
	b = unlist(post["b[1]"], use.names=F),
	stoll_factor= unlist(post["ec_factor"], use.names=F),
	expulsion_state = e.vec,
	autopsy_state = a.vec,
	log_lik=log_lik)

#Export as an r object
saveRDS(posterior.list, file = "posteriors/wm-pl-nbh-2p.rds")
