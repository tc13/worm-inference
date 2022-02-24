//inference model for worm burden from expulsion studies
//Multi-study model, power law function
//Negative binomial egg error

data {
  int<lower=1> N_expul;     //Number of individuals in worm expulsion studies studies
  int<lower=1> N_autopsy;  //Number of individuals in autopsy studies study
  int<lower=0> epg_expul[N_expul];     //Reported eggs per gram per person (expulsion)
  int<lower=0> worms_expul[N_expul];   //Observed worms expelled per person (expulsion)
  int<lower=0> epg_autopsy[N_autopsy];     //Reported eggs per gram per person (autopsy)
  int<lower=0> worms_autopsy[N_autopsy];   //Observed worms expelled per person (autopsy)
  int<lower=2> N_expul_studies;  //Number of expulsion studies
  int<lower=1> N_autopsy_studies; //Num of autopsy studies
  int<lower=1, upper=N_expul_studies> study_id[N_expul]; //study ID
  int<lower=1,upper=N_expul_studies> ramsay_id; //indicates which study is Ramsay et al. using Stoll dilution method
  int<lower=1> delta_worm; //Difference between observed worm burden and possible max
}

transformed data{
  int<lower=1> N;
  int<lower=2> N_studies;
  int max_worm;
  N = (N_expul+N_autopsy);
  N_studies = (N_expul_studies+N_autopsy_studies);
  max_worm = max(worms_autopsy)+delta_worm;
}

parameters {
  real<lower=0> y1; //worm fecundity param (power law)
  real<lower=0,upper=1> gamma; //worm fecundity param (power law)
  real<lower=0> M[N_studies];   //Mean worm burden
  real<lower=0> k[N_studies];  //dispersion of worms
  real<lower=0, upper=1> pr_recovery; //probability of worm recovery from expulsion
  real<lower=0> k_mean; //hyper-parameter for k
  real<lower=0> k_sd; //hyper-parameter for k
  real<lower=1> stoll_factor;
  real<lower=0> h; //neg binom variance in egg output
}

transformed parameters{
   matrix[N_expul,delta_worm] marginal_expul; //Marginal probability of each possible worm value
   matrix[N_autopsy, 4] marginal_autopsy; //Marginal probability of each possible worm value
   vector[max_worm] epg_expected; //Expected egg output given worm value
   row_vector[N] epg_factor; //EPG correction factor
   epg_factor = rep_row_vector(1, N); //defaults to 1
   for(i in 1:max_worm){
     epg_expected[i] = y1*i^gamma; //power law function;
   }
   //Expulsion studies likelihood
   for(i in 1:N_expul){
     if(study_id[i]==ramsay_id)
        epg_factor[i] = stoll_factor; //multiplication factor for Ramsay study
     
     if(epg_expul[i]==0 && worms_expul[i]==0){ //if no observed eggs or worms
          marginal_expul[i,1] = neg_binomial_2_lpmf(0 | M[study_id[i]], k[study_id[i]]);
          for(j in 2:delta_worm){
              marginal_expul[i,j] = neg_binomial_2_lpmf(0 | (epg_expected[(worms_expul[i]+j-1)]/epg_factor[i]),h) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery);
          }
      }
        else if(epg_expul[i]>0 && worms_expul[i]==0){ //if eggs observed but no worms
           marginal_expul[i,1] = negative_infinity(); //impossible that there are zero worms
           for(j in 2:delta_worm){
              marginal_expul[i,j] = neg_binomial_2_lpmf(epg_expul[i] | (epg_expected[(worms_expul[i]+j-1)]/epg_factor[i]), h) + neg_binomial_2_lpmf((j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(0 | (j-1), pr_recovery); 
           }
        }else{ //if >0 worms observed
          for(j in 1:delta_worm){
            marginal_expul[i,j] = neg_binomial_2_lpmf(epg_expul[i] | (epg_expected[(worms_expul[i]+j-1)]/epg_factor[i]), h) + neg_binomial_2_lpmf((worms_expul[i]+j-1) | M[study_id[i]], k[study_id[i]]) + binomial_lpmf(worms_expul[i] | (worms_expul[i]+j-1), pr_recovery);
          }
        }
      }
    //Autopsy study likelihood - can divide epg by factor of 1,2,3 or 4
    for(i in 1:N_autopsy){
        if(epg_autopsy[i]==0 && worms_autopsy[i]==0)
          marginal_autopsy[i] = rep_row_vector(log(0.25),4);
        else{
          for(j in 1:4){
            marginal_autopsy[i,j] = neg_binomial_2_lpmf(epg_autopsy[i] | epg_expected[worms_autopsy[i]]*j, h) + neg_binomial_2_lpmf(worms_autopsy[i] | M[N_studies], k[N_studies]);
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
    y1 ~ gamma(20, 5);
    gamma ~ beta(200, 200); //strong prior for density dependence - informed by extracted worm fecund. data
    M[1] ~ normal(39, 20);  //prior for Elkins study
    M[2] ~ normal(187, 20); //prior for Sayasone study
    M[3] ~ normal(85, 20);  //prior for Ramsay study
    M[4] ~ normal(160, 20); //prior for Autopsy study
    k ~ normal(k_mean, k_sd);
    pr_recovery ~ beta(850,200);
    k_mean ~ normal(0.5, 2);
    k_sd ~ normal(0.5, 1);
    h ~ normal(10, 2);
    stoll_factor ~ normal(100, 50);
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
