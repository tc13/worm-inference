// Inference of (liver fluke) worm burden from egg counts
// Tom Crellen, June 2021, thomas.crellen@bdi.ox.ac.uk 
// Adapted from de Vlas et al. 1992 https://doi.org/10.1017/S003118200006371X
// Multiple egg counts per-person
// Worm burden fitted on log scale
// Incorporates Antigen test results
// Multiple clusters - variable M by cluster - k is constant
// Negative binomial, density dependent egg output (power-law)

data {
  int<lower=1> N; //individuals
  int<lower=1> K; //total egg count observations
  int<lower=1> N_clusters; //number of clusters (random effects)
  int<lower=1, upper=N_clusters> cluster[N]; //cluster of each individual
  int<lower=1> slides[N]; //observations per person
  int<lower=1> slide1[N]; //index first slide per person
  int<lower=1> slide2[N]; //index second slide per person
  int<lower=0> eggs[K]; //egg counts
  int<lower=0, upper=1> urine[N]; //indicates if urine antigen test taken
  int<lower=0, upper=1> antigen[N]; //antigen test result
  int<lower=100> big_int; //sum values of n from zero to this number
  int<lower=0, upper=1> egg_neg[N]; //indicates if person egg negative or not
  real<lower=0> stool_mass[N]; //mass of stool sample in grams
   int<lower=1> stool_drops[N]; //number of drops in FECT concentrate
  //worm fecundity data
  int<lower=0> Nw;
  int<lower=0> E[Nw];
  int<lower=0> W[Nw];
  int<lower=0> N_studies;
  int<lower=0, upper=N_studies> study[Nw];
}

transformed data{
  row_vector[big_int] worms;
  row_vector[N] stool_factor;
  int reffects;
  for(j in 1:big_int)
    worms[j] = j-1;
  //factor to modify EPG to eggs
  for(i in 1:N)
    stool_factor[i] = stool_mass[i]/stool_drops[i];
  reffects = N_clusters+N_studies;
}  

parameters {
  real M_ln[reffects]; //mean log worm burden in population
  real<lower=0> k; //worm dispersion parameter
  real<lower=0, upper=1> sp; //antigen test specificity
  real<lower=0, upper=1> se; //antigen test sensitivity
  real M_mu; //Log normal hyperprior of M
  real<lower=0> M_sd;  //Log normal hyperprior of M
  real<lower=0> y1; //worm fecundity param (power law)
  real<lower=0,upper=1> gamma; //worm fecundity param (power law)
}

transformed parameters {
  real<lower=0> M[reffects];
  real<lower=0, upper=1> p[reffects]; //prevalence from negative binomial distribution
  M = exp(M_ln);
  for(i in 1:reffects)
    p[i] = 1-(1+M[i]/k)^-k; //calculate p for each cluster
}

model {
  int counter; //track egg index
  row_vector[big_int] epg;  //store epg for integer worm counts
  row_vector[N_clusters] antigen_prob;
  matrix[N_clusters, big_int] worm_prior;
  vector[Nw] expected_epg;  //store epg predicted by model
  
  //worm fecundity data
  for(i in 1:Nw){
    expected_epg[i] = y1*W[i]^gamma; //power law function
    if (E[i]==0 && W[i]==0)
        target += 1;
    else if (E[i]>0 && W[i]==0)
        target += negative_infinity();
    else
        target += poisson_lpmf(E[i] | expected_epg[i]) + neg_binomial_2_lpmf(W[i] | M[(study[i]+N_clusters)], k); 
  }
  
  //predicted epg given worm burden (starts at zero)
  for(i in 1:big_int)
    epg[i] = y1*worms[i]^gamma;
  
  //Calculate prob. antigen test positive per cluster
  for(i in 1:N_clusters){
    antigen_prob[i] = (1-sp) + p[i]*(se+sp-1);
    //Compute worm prior per cluster
    for(j in 1:big_int){
      worm_prior[i,j] = neg_binomial_2_lpmf((j-1) | M[i], k);
    }
  }
  
//Calculate probabilites
  for(i in 1:N){
    row_vector[big_int] marginal;
    row_vector[big_int] expected_eggs;
    expected_eggs = epg*stool_factor[i]; //eggs per gram * grams of stool / number of drops
    if(egg_neg[i]==0){ //if all egg counts are zero
      marginal[1] = worm_prior[cluster[i],1] + bernoulli_lpmf(antigen[i] | antigen_prob[cluster[i]])*urine[i]; //special case, worms (j)=0
      for(j in 2:big_int){  //worms (j) from 1:(big_int-1)
            marginal[j] = poisson_lpmf(0 | expected_eggs[j])*slides[i] + worm_prior[cluster[i], j] + bernoulli_lpmf(antigen[i] | antigen_prob[cluster[i]])*urine[i];
       }
     }else{ //at least 1 non-zero egg count
      marginal[1] =  negative_infinity(); //special case, worms (j)=0
      for(j in 2:big_int){  //worms (j) from 1:(big_int-1)
            marginal[j] = poisson_lpmf(eggs[slide1[i]:slide2[i]] | expected_eggs[j]) + worm_prior[cluster[i], j] + bernoulli_lpmf(antigen[i] | antigen_prob[cluster[i]])*urine[i];
        }
     }

   //increment log likelihood
   target += log_sum_exp(marginal);
  }
  
  //prior distributions
  y1 ~ gamma(20, 5);
  gamma ~ beta(755, 245); //strong prior for density dependence - informed by extracted worm fecund. data
  M_ln ~ normal(M_mu, M_sd);
  M_mu ~ normal(0, 4);    //allow M to vary by cluster
  M_sd ~ normal(1, 3);
  k ~ normal(0.1, 0.2);   //strong prior on k - Burli et al.
  sp ~ beta(669, 65);     //Strong priors - from Worasith et al. 2019 https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0007186 
  se ~ beta(269, 44); 
}
