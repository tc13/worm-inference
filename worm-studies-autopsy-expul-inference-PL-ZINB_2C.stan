//Inference model for worm burden from expulsion studies & autopsy
//Multi-study model, power-law function
//Zero inflated neg binom variance in egg output
//Two functions

data {
  int<lower=1> N_expul;     //Number of individuals in worm expulsion studies studies
  int<lower=1> N_autopsy;  //Number of individuals in autopsy studies study
  int<lower=0> epg_expul[N_expul];     //Reported eggs per gram per person (expulsion)
  int<lower=0> worms_expul[N_expul];   //Observed worms expelled per person (expulsion)
  int<lower=0> epg_autopsy[N_autopsy];     //Reported eggs per gram per person (autopsy)
  int<lower=0> worms_autopsy[N_autopsy];   //Observed worms expelled per person (autopsy)
  int<lower=2> N_expul_studies;  //Number of expulsion studies
  int<lower=1, upper=N_expul_studies> study_id[N_expul]; //study ID
  int<lower=1> delta_worm; //Difference between observed worm burden and possible max
  int<lower=0, upper=N_expul_studies> stoll_study; //which study ID uses Stoll diagnostic, if none can be set to zero
  int<lower=0, upper=N_expul_studies> lao_study; //which study is from Lao (uses seperate relationship) 
}

transformed data{
  int<lower=1> N;
  int<lower=2> N_studies;
  int max_worm;
  N = N_expul+N_autopsy;
  N_studies = N_expul_studies+1;
  max_worm = max(worms_autopsy)+delta_worm;
}

parameters {
  real<lower=0> y1[2]; //worm fecundity param (power law)
  real<lower=0,upper=1> gamma[2]; //worm fecundity param (power law)
  real<lower=0> h[2]; //neg binom variance in egg output
  real<lower=0> h_mean; //hyper parameter for h
  real<lower=0> h_sd; //hyper parameter for h
  real<lower=0> M[N_studies];   //Mean worm burden
  real<lower=0> k[N_studies];  //dispersion of worms
  real<lower=0, upper=1> pr_recovery; //probability of worm recovery from expulsion
  real<lower=0> k_mean; //hyper-parameter for k
  real<lower=0> k_sd; //hyper-parameter for k
  real<lower=0> sens_b; //Michaelis-Menten sensitivity parameter (ZINB)
  real<lower=1> stoll_factor; //factor to multiply Stoll by
}

transformed parameters{
   matrix[N_expul,delta_worm] marginal_expul; //Marginal probability of each possible worm value
   matrix[N_autopsy, 4] marginal_autopsy; //Marginal probability of each possible worm value
   vector[max_worm] epg_expected_TH; //Expected egg output given worm value (THAILAND)
   vector[max_worm] epg_expected_LAO; //Expected egg output given worm value (LAOS)
   real sens[max_worm]; //sensitivity given worm burden (catalytic function)
   for(i in 1:max_worm){
     epg_expected_TH[i] = y1[1]*i^gamma[1]; //power law function (THAILAND)
     epg_expected_LAO[i] = y1[2]*i^gamma[2]; //power law function (LAOS)
     sens[i] = i/(sens_b+i);   
     }
   //Expulsion studies likelihood
   for(i in 1:N_expul){
      if(epg_expul[i]==0 && worms_expul[i]==0){ //if no eggs or worms observed
        marginal_expul[i,1] = neg_binomial_2_lpmf(0 | M[study_id[i]], k[study_id[i]]); //true negative worms and eggs
        for(j in 2:delta_worm){ //false negative worms and eggs
            marginal_expul[i,j] = bernoulli_lpmf(0 | sens[(j-1)]) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery);
          }
      }
      else if(epg_expul[i]==0 && worms_expul[i]>0){ //worms but no eggs (false negative eggs)
          for(j in 1:delta_worm){
            marginal_expul[i,j] = bernoulli_lpmf(0 | sens[(worms_expul[i]+j-1)]) + neg_binomial_2_lpmf((worms_expul[i]+j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(worms_expul[i] | (worms_expul[i]+j-1), pr_recovery);
          }
        }
        else if(epg_expul[i]>0 && worms_expul[i]==0){  //if eggs observed but no worms
           marginal_expul[i,1] = negative_infinity();  //impossible that there are zero worms
           
          if(stoll_study==study_id[i]){
            for(j in 2:delta_worm){
              marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | (epg_expected_TH[(j-1)]*stoll_factor), h[1]) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery);
            }
          }else if(lao_study==study_id[i]){
           for(j in 2:delta_worm){
              marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | epg_expected_LAO[(j-1)], h[2]) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery); 
           }
          }else{
              for(j in 2:delta_worm){
                marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | epg_expected_TH[(j-1)], h[1]) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery); 
           }
        }
        }
        else{  //if >0 worms and eggs observed
            if(stoll_study==study_id[i]){
               for(j in 1:delta_worm){
                  marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(worms_expul[i]+j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | (epg_expected_TH[(worms_expul[i]+j-1)]*stoll_factor), h[1]) + neg_binomial_2_lpmf((worms_expul[i]+j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(worms_expul[i] | (worms_expul[i]+j-1), pr_recovery);
                }
             }else if(lao_study==study_id[i]){
               for(j in 1:delta_worm){
                   marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(worms_expul[i]+j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | epg_expected_LAO[(worms_expul[i]+j-1)], h[2]) + neg_binomial_2_lpmf((worms_expul[i]+j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(worms_expul[i] | (worms_expul[i]+j-1), pr_recovery);
               }
             }else{
                for(j in 1:delta_worm){
                   marginal_expul[i,j] = bernoulli_lpmf(1 | sens[(worms_expul[i]+j-1)]) + neg_binomial_2_lpmf(epg_expul[i] | epg_expected_TH[(worms_expul[i]+j-1)], h[1]) + neg_binomial_2_lpmf((worms_expul[i]+j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(worms_expul[i] | (worms_expul[i]+j-1), pr_recovery);
                }
              }
          }
      }
    //Autopsy study likelihood - can divide epg by factor of 1,2,3 or 4
    for(i in 1:N_autopsy){
        if(epg_autopsy[i]==0 && worms_autopsy[i]==0) //no eggs or worms
          marginal_autopsy[i] = rep_row_vector(log(0.25),4);
        else if(epg_autopsy[i]==0 && worms_autopsy[i]>0){ //eggs false negative
          for(j in 1:4){
            marginal_autopsy[i,j] = bernoulli_lpmf(0 | sens[worms_autopsy[i]]) + neg_binomial_2_lpmf(worms_autopsy[i] | M[N_studies], k[N_studies]);
          }
        }else{ //positive for eggs and worms
          for(j in 1:4){
            marginal_autopsy[i,j] = bernoulli_lpmf(1 | sens[worms_autopsy[i]]) + neg_binomial_2_lpmf(epg_autopsy[i] | epg_expected_TH[worms_autopsy[i]]*j, h[1]) + neg_binomial_2_lpmf(worms_autopsy[i] | M[N_studies], k[N_studies]);
          }
        }
    }
}

model{
    //increment log likelihood for expulsion studies
    for(i in 1:N_expul)
      target += log_sum_exp(marginal_expul[i]);
    //increment log likelihood for autopsy study
    for(i in 1:N_autopsy)
      target += log_sum_exp(marginal_autopsy[i]);
    
    //prior distributions
    y1 ~ normal(12, 2);     //expected output from 1 worm: 3160 (Wykoff & Ariyaprakai, Opisthorchis viverrini in Thailand-egg production in man and laboratory animals. Journal of Parasitology 52:4 (1966)) divided by daily mass of human stool - 250g for developing countries (Rose, C., Parker, A., Jefferson, B. and Cartmell, E., 2015. The characterization of feces and urine: a review of the literature to inform advanced treatment technology. Critical reviews in environmental science and technology, 45(17), pp.1827-1879)
    gamma ~ beta(20, 20);   //prior for density dependence
    M[1] ~ normal(39, 10);  //prior for Elkins study
    M[2] ~ normal(187, 10); //prior for Sayasone study
    M[3] ~ normal(85, 10);  //prior for Ramsay study
    M[4] ~ normal(160, 10); //prior for Autopsy study
    M[5] ~ normal(49, 10); //prior for Haswell study
    k ~ normal(k_mean, k_sd); //heriarchical values for k
    k_mean ~ exponential(1);  //mean of k
    k_sd ~ exponential(2); //variance of k 
    pr_recovery ~ beta(2, 2); //probability of worm recovery
    h ~ normal(h_mean, h_sd);  //egg count dispersion
    h_mean ~ exponential(0.2); //mean of h
    h_sd ~ exponential(2); //sd of h
    sens_b ~ normal(1, 1); //egg count sensitivity parameter
    stoll_factor ~ normal(50, 5); //stoll factor parameter
}

generated quantities{
  matrix[N_expul,delta_worm] pState; //normalised probabilties for latent discrete parameter
  matrix[N_autopsy, 4] aut_state;    //normalised probabilties for latent discrete parameter
  vector[N] log_lik;    //store log likelihood for model comparison 
  for(i in 1:N_expul){
    pState[i] = exp(marginal_expul[i] - log_sum_exp(marginal_expul[i])); //normalise marginal prob
    log_lik[i] = log_sum_exp(marginal_expul[i]);
  }
  for(i in 1:N_autopsy){
    aut_state[i] = exp(marginal_autopsy[i] - log_sum_exp(marginal_autopsy[i]));
    log_lik[(i+N_expul)] = log_sum_exp(marginal_autopsy[i]);
  }
}
