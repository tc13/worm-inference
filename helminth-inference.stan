//Inference model for worm burden from expulsion & autopsy studies
//Flexible for egg output function and error distribution
//Multiple populations with seperate density dependent functions
//Thomas Crellen, University of Glasgow, September 2022, thomas.crellen@glasgow.ac.uk 
data {
  int<lower=1> N; //total participants
  array[N] int<lower=0> epg; //reported eggs counts per person
  array[N] int<lower=0> worms; //recovered worms per person
  int<lower=1> N_studies; //total number of studies
  array[N] int<lower=1, upper=N_studies> study_id; //index individual to study
  array[N_studies] real<lower=0> mean_worms_obs; //mean number of worms recovered per study (for prior)
  int<lower=1> N_pops; //Number of populations with seperate worm/egg relationship
  array[N_studies] int<lower=1> study_pop; //match study to populations with seperate worm/egg relationship
  int<lower=1, upper=N> N_expul; //number of individuals in expulsion studies
  int<lower=1, upper=N> N_autopsy; //number of individuals in autopsy studies
  array[N_expul] int<lower=1, upper=N> expul_indx; //indices for individuals in expulsion study
  array[N_autopsy] int<lower=1, upper=N> autopsy_indx; //indices for individuals in expulsion study
  array[N_studies] int<lower=0, upper=1> ecc_id; //variable if study needs egg count correction
  int<lower=1, upper=2> model_function; //1 = Power Law function | 2 = Algebraic Decay function
  int<lower=1, upper=3> error_dist; //1 = Poisson | 2 = Negative Binomial | 3 = Negative Binomial Hurdle
  int<lower=1> delta_worm; //difference between observed worm burden and possible max
}
transformed data {
  array[N_expul] int<lower=0> epg_expul; //Reported eggs per gram per person (expulsion)
  array[N_expul] int<lower=0> worms_expul; //Observed worms expelled per person (expulsion)
  array[N_autopsy] int<lower=0> epg_autopsy; //Reported eggs per gram per person (autopsy)
  array[N_autopsy] int<lower=0> worms_autopsy; //Observed worms expelled per person (autopsy)
  array[N_expul] int<lower=1, upper=N_studies> expul_study_id;
  array[N_autopsy] int<lower=1, upper=N_studies> autopsy_study_id;
  int max_worm;
  array[N_expul] int<lower=1, upper=N_pops> pop_id_expul;
  array[N_autopsy] int<lower=1, upper=N_pops> pop_id_autopsy;
  array[N_expul] int<lower=1, upper=4> cases_expul;
  array[N_autopsy] int<lower=1, upper=4> cases_autopsy;
  epg_expul = epg[expul_indx];
  worms_expul = worms[expul_indx];
  epg_autopsy = epg[autopsy_indx];
  worms_autopsy = worms[autopsy_indx];
  expul_study_id = study_id[expul_indx];
  autopsy_study_id = study_id[autopsy_indx];
  max_worm = max(worms_expul) + delta_worm;
  //Set case conditions
  for (i in 1 : N_expul) {
    pop_id_expul[i] = study_pop[expul_study_id[i]];
    if (epg_expul[i] == 0 && worms_expul[i] == 0)  //no eggs or worms observed (true negative)
      cases_expul[i] = 1;
    else if (epg_expul[i] > 0 && worms_expul[i] == 0)  //eggs observed but no worms (false negative worms)
      cases_expul[i] = 2;
    else if (epg_expul[i] == 0 && worms_expul[i] > 0)  //worms observed but no eggs (false negative eggs)
      cases_expul[i] = 3;
    else 
      cases_expul[i] = 4; //eggs and worms detected (true positive)
  }
  for (i in 1 : N_autopsy) {
    pop_id_autopsy[i] = study_pop[autopsy_study_id[i]];
    if (epg_autopsy[i] == 0 && worms_autopsy[i] == 0)  //no eggs or worms observed (true negative)
      cases_autopsy[i] = 1;
    else if(epg_autopsy[i] == 0 && worms_autopsy[i] > 0)
      cases_autopsy[i] = 2; //worms observed but no eggs (false negative eggs)
    else
      cases_autopsy[i] = 3; //eggs and worms detected (true positive)
  }
}
parameters {
  array[N_studies] real<lower=0> M; //Mean worm burden
  array[N_studies] real<lower=0> k; //dispersion of worms
  real<lower=0> k_mean; //hyper-parameter for k
  real<lower=0> k_sd; //hyper-parameter for k
  real<lower=0, upper=1> pr_recovery; //probability of worm recovery from expulsion
  real<lower=1> ec_factor; //factor to multiply mean egg counts by
  array[(model_function == 1) ? N_pops : 0] real<lower=0> y1; //worm fecundity param (power law)
  array[(model_function == 1) ? N_pops : 0] real<lower=0, upper=1> gamma; //worm fecundity param (power law)
  array[(model_function == 2) ? N_pops : 0] real<lower=0> L0; //worm fecundity param (algebraic decay)
  array[(model_function == 2) ? N_pops : 0] real<lower=0> M0; //worm fecundity param (algebraic decay)
  array[(error_dist > 1) ? N_pops : 0] real<lower=0> h; //neg binom variance in egg output
  array[(error_dist == 3) ? 1 : 0] real<lower=0> b; //Michaelis-Menten sensitivity parameter (NBH)
}
transformed parameters {
  matrix[N_expul, delta_worm] marginal_expul; //Marginal probability of each possible worm value
  matrix[N_autopsy, 4] marginal_autopsy; //Marginal probability of each possible worm value
  matrix[max_worm, N_studies] epg_expected; //Expected egg output given worm value x
  array[(error_dist == 3) ? max_worm : 0] real sens; //sensitivity given worm burden (catalytic function)
  vector[N_studies] ec_study; //vector for correcting egg counts into epg by study
  
  //Vector for egg count correction by population
  ec_study = rep_vector(1, N_studies);
  for (j in 1 : N_studies) {
    if (ecc_id[j] == 1) {
      ec_study[j] = ec_factor;
    }
  }
  
  //Expected EPG depends on function
  if (model_function == 1) {
    for (x in 1 : max_worm) {
      for (j in 1 : N_studies) {
        epg_expected[x, j] = (y1[study_pop[j]]*x^gamma[study_pop[j]] ) * ec_study[j]; //power law function
      }
    }
  } 
  
  if (model_function == 2) {
    for (x in 1 : max_worm) {
      for (j in 1 : N_studies) {
        epg_expected[x, j] = ((x * L0[study_pop[j]] * M0[study_pop[j]]) / (x + M0[study_pop[j]])) * ec_study[j]; //algebraic decay function
      }
    }
  }
  //Hurdle component for Error Dist 3 (NBH)
  if (error_dist == 3) {
    for (x in 1 : max_worm) 
      sens[x] = x / (b[1] + x);
  }
  
  //Expulsion studies likelihood
  //Error distribution 1 = Poisson
  if (error_dist == 1) {
    for (i in 1 : N_expul) {
      if (cases_expul[i] == 1) {
        marginal_expul[i, 1] = neg_binomial_2_lpmf(0 | M[study_id[i]], k[study_id[i]]); //true negative worms and eggs
        for (j in 2 : delta_worm) { //false negative worms and eggs
          marginal_expul[i, j] = poisson_lpmf(0 | (epg_expected[(j - 1), expul_study_id[i]]))
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else if (cases_expul[i] == 2) {
        //if eggs observed but no worms
        marginal_expul[i, 1] = negative_infinity(); //impossible that there are zero worms
        for (j in 2 : delta_worm) {
          marginal_expul[i, j] = poisson_lpmf(epg_expul[i] | (epg_expected[(j - 1), expul_study_id[i]]))
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else {
        //if >0 worms observed
        for (j in 1 : delta_worm) {
          marginal_expul[i, j] = poisson_lpmf(epg_expul[i] | (epg_expected[(worms_expul[i]  + j - 1), expul_study_id[i]]))
                                 + neg_binomial_2_lpmf((worms_expul[i] + j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(worms_expul[i] | (worms_expul[i] + j - 1), pr_recovery);
        }
      }
    }
    //Error distribution 2 = Negative Binomial
  } else if (error_dist == 2) {
    for (i in 1 : N_expul) {
      if (cases_expul[i] == 1) {
        marginal_expul[i, 1] = neg_binomial_2_lpmf(0 | M[expul_study_id[i]], k[expul_study_id[i]]);
        for (j in 2 : delta_worm) {
          marginal_expul[i, j] = neg_binomial_2_lpmf(0 | (epg_expected[(j - 1), expul_study_id[i]]), h[pop_id_expul[i]])
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else if (cases_expul[i] == 2) {
        //if eggs observed but no worms
        marginal_expul[i, 1] = negative_infinity(); //impossible that there are zero worms
        for (j in 2 : delta_worm) {
          marginal_expul[i, j] = neg_binomial_2_lpmf(epg_expul[i] | (epg_expected[(j - 1), expul_study_id[i]]), h[pop_id_expul[i]])
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else {
        //if >0 worms observed
        for (j in 1 : delta_worm) {
          marginal_expul[i, j] = neg_binomial_2_lpmf(epg_expul[i] | (epg_expected[(worms_expul[i] + j - 1), expul_study_id[i]]), h[pop_id_expul[i]])
                                 + neg_binomial_2_lpmf((worms_expul[i] + j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(worms_expul[i] | (worms_expul[i] + j - 1), pr_recovery);
        }
      }
    }
    //Error distribution 3 = Negative Binomial Hurdle
  } else {
    for (i in 1 : N_expul) {
      if (cases_expul[i] == 1) {
        marginal_expul[i, 1] = neg_binomial_2_lpmf(0 | M[expul_study_id[i]], k[expul_study_id[i]]);
        for (j in 2 : delta_worm) { //false negative worms and eggs
          marginal_expul[i, j] = bernoulli_lpmf(0 | sens[(j - 1)])
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else if (cases_expul[i] == 2) {
        marginal_expul[i, 1] = negative_infinity(); //impossible that there are zero worms
        for (j in 2 : delta_worm) {
          marginal_expul[i, j] = bernoulli_lpmf(1 | sens[(j - 1)])
                                 + neg_binomial_2_lpmf(epg_expul[i] | epg_expected[(j - 1), expul_study_id[i]], h[pop_id_expul[i]])
                                 - neg_binomial_2_lccdf(0 | epg_expected[(j - 1), expul_study_id[i]], h[pop_id_expul[i]])
                                 + neg_binomial_2_lpmf((j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(0 | (j - 1), pr_recovery);
        }
      } else if (cases_expul[i] == 3) {            //observe worms but no eggs
        for (j in 1 : delta_worm) {
          marginal_expul[i, j] = bernoulli_lpmf(0 | sens[(worms_expul[i] + j - 1)])
                                 + neg_binomial_2_lpmf((worms_expul[i] + j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(worms_expul[i] | (worms_expul[i] + j - 1), pr_recovery);
        }
      } else {
        for (j in 1 : delta_worm) {
          marginal_expul[i, j] = bernoulli_lpmf(1 | sens[(worms_expul[i] + j - 1)])
                                 + neg_binomial_2_lpmf(epg_expul[i] | epg_expected[(worms_expul[i] + j - 1), expul_study_id[i]], h[pop_id_expul[i]])
                                 - neg_binomial_2_lccdf(0 | epg_expected[(worms_expul[i] + j - 1), expul_study_id[i]], h[pop_id_expul[i]])
                                 + neg_binomial_2_lpmf((worms_expul[i] + j - 1) | M[expul_study_id[i]], k[expul_study_id[i]])
                                 + binomial_lpmf(worms_expul[i] | (worms_expul[i] + j - 1), pr_recovery);
        }
      }
    }
  }
  //Autopsy study likelihood - divide epg by factor of 1,2,3 or 4
  //Error distribution 1 = Poisson
  if (error_dist == 1) {
    for (i in 1 : N_autopsy) {
      if (cases_autopsy[i] == 1) {  //no eggs or worms
        marginal_autopsy[i] = rep_row_vector(log(0.25), 4);
      }
      else {
        //worms recovered
        for (j in 1 : 4) {
          marginal_autopsy[i, j] = poisson_lpmf(epg_autopsy[i] | epg_expected[worms_autopsy[i], autopsy_study_id[i]]*j)
                                   + neg_binomial_2_lpmf(worms_autopsy[i] | M[autopsy_study_id[i]], k[autopsy_study_id[i]]);
        }
      }
    }
    //Error distribution 2 = Negative Binomial
  } else if (error_dist == 2) {
    for (i in 1 : N_autopsy) {
      if (cases_autopsy[i] == 1) { //no eggs or worms
        marginal_autopsy[i] = rep_row_vector(log(0.25), 4);
      }
      else {
        //worms recovered
        for (j in 1 : 4) {
          marginal_autopsy[i, j] = neg_binomial_2_lpmf(epg_autopsy[i] | epg_expected[worms_autopsy[i], autopsy_study_id[i]]*j, h[pop_id_autopsy[i]])
                                   + neg_binomial_2_lpmf(worms_autopsy[i] | M[autopsy_study_id[i]], k[autopsy_study_id[i]]);
        }
      }
    }
  } else {
    //Error distribution 3 = Negative Binomial Hurdle
    for (i in 1 : N_autopsy) {
      if (cases_autopsy[i] == 1) { //no eggs or worms
        marginal_autopsy[i] = rep_row_vector(log(0.25), 4);
      }
      else if (cases_autopsy[i] == 2) {
        //eggs false negative
        for (j in 1 : 4) {
          marginal_autopsy[i, j] = bernoulli_lpmf(0 | sens[worms_autopsy[i]])
                                   + neg_binomial_2_lpmf(worms_autopsy[i] | M[autopsy_study_id[i]], k[autopsy_study_id[i]]);
        }
      } else {
        //positive for eggs and worms
          for (j in 1 : 4) {
           marginal_autopsy[i, j] = bernoulli_lpmf(1 | sens[worms_autopsy[i]])
                                   + neg_binomial_2_lpmf(epg_autopsy[i] | epg_expected[worms_autopsy[i], autopsy_study_id[i]]*j, h[pop_id_autopsy[i]])
                                   - neg_binomial_2_lccdf(0 | epg_expected[worms_autopsy[i], autopsy_study_id[i]]*j, h[pop_id_autopsy[i]])
                                   + neg_binomial_2_lpmf(worms_autopsy[i] | M[autopsy_study_id[i]], k[autopsy_study_id[i]]);
        }
      }
    }
  }
}
model {
  //increment log likelihood for expulsion studies
  for (i in 1 : N_expul) {
    target += log_sum_exp(marginal_expul[i]);
  }
  //increment log likelihood for autopsy study
  for (i in 1 : N_autopsy) {
    target += log_sum_exp(marginal_autopsy[i]);
  }
  
  //prior distributions
  y1 ~ normal(20, 5); //expected output from 1 worm: 3160 (Wykoff & Ariyaprakai, Opisthorchis viverrini in Thailand-egg production in man and laboratory animals. Journal of Parasitology 52:4 (1966)) divided by daily mass of human stool - 250g for developing countries (Rose, C., Parker, A., Jefferson, B. and Cartmell, E., 2015. The characterization of feces and urine: a review of the literature to inform advanced treatment technology. Critical reviews in environmental science and technology, 45(17), pp.1827-1879)
  gamma ~ beta(4, 4); //prior for density dependence
  L0 ~ normal(20, 5); //expected output from 1 worm: 3160 (Wykoff & Ariyaprakai, Opisthorchis viverrini in Thailand-egg production in man and laboratory animals. Journal of Parasitology 52:4 (1966)) divided by daily mass of human stool - 250g for developing countries (Rose, C., Parker, A., Jefferson, B. and Cartmell, E., 2015. The characterization of feces and urine: a review of the literature to inform advanced treatment technology. Critical reviews in environmental science and technology, 45(17), pp.1827-1879)
  M0 ~ normal(2000, 100);
  
  for (j in 1 : N_studies) {
    M[j] ~ normal(mean_worms_obs[j], 20); //recovered mean worm burden per study
  }
  k ~ normal(k_mean, k_sd);
  k_mean ~ exponential(1);
  k_sd ~ exponential(5);
  pr_recovery ~ beta(4, 2);
  h ~ gamma(10, 4);
  b ~ normal(2, 0.1); //egg count sensitivity parameter - informed by Lovis et al.
  ec_factor ~ normal(5, 10); //stoll factor parameter
}
generated quantities {
  matrix[N_expul, delta_worm] pState; //normalised probabilties for latent discrete parameter
  matrix[N_autopsy, 4] aut_state; //normalised probabilties for latent discrete parameter
  vector[N] log_lik; //store log likelihood for model comparison 
  for (i in 1 : N_expul) {
    pState[i] = exp(marginal_expul[i] - log_sum_exp(marginal_expul[i])); //normalise marginal prob
    log_lik[i] = log_sum_exp(marginal_expul[i]);
  }
  for (i in 1 : N_autopsy) {
    aut_state[i] = exp(marginal_autopsy[i] - log_sum_exp(marginal_autopsy[i]));
    log_lik[(i + N_expul)] = log_sum_exp(marginal_autopsy[i]);
  }
}
