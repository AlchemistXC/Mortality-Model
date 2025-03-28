functions {
  real partial_sum_lpdf(array[] real logm,
                        int start, int end,
                        array[] int death) {
    return poisson_log_lupmf(death[start:end] | logm);
  }
}
data {
  int<lower=0> T;                 // number of years
  int<lower=0> A;                  // number of age categories
  int<lower=0> P;                  // number of age categories
  int<lower=0> Tf;                 // number of forecast years
  int<lower=0> AT;
  int<lower=0> ATP;
  int<lower=0> ATfP;
  array[ATP] real loge;        // log exposures
  array[ATP] int<lower=0> death;   // deaths
  vector[A*P] mu_a;
}
transformed data{
  int grainsize = 1;
}
parameters {
  vector[A*P] ap;                    // alpha_x
  array[P] simplex[A] bp; 
  array[P] real cp;
  array[P] vector[T-1] kps;                 // vector of k_t differences
  array[P] real<lower = 0> sigma_kp;          // standard deviation of the random walk
}
transformed parameters  {
  array[P] vector[T] kp; 
  array[ATP] real logm;
  
  for(p in 1:P){
    kp[p][1] = 0;
    //kp[p][2] = -sum(kps[p][1:(T-2)]);
    kp[p][2:T] = kps[p];
  }
  for(p in 1:P){
    for (j in 1:T){
      for (i in 1:A) {
      logm[AT*(p-1) + A * (j-1) + i] = loge[AT*(p-1) + A * (j-1) + i] + ap[(p-1)*A + i] + bp[p][i] * kp[p][j];
      }
    }
  }
  
  
}
model {

  target += reduce_sum(partial_sum_lpdf, logm, grainsize, death);
  //target += poisson_log_lpmf(death| logm);    
  for(p in 1:P){
    target += normal_lpdf(kps[p][1]| cp[p], sigma_kp[p]);
    target += normal_lpdf(kps[p][2:(T-1)]| cp[p] + kps[p][1:(T-2)], sigma_kp[p]);
    target += dirichlet_lpdf(bp[p]|rep_vector(1, A));
    target += normal_lpdf(cp[p]|0, sqrt(10));
    //target += normal_lpdf(sigma_k|2, 1);
    //target += student_t_lupdf(sigma_kp[p]|4, 0, 1);
    //target += exponential_lpdf(sigma_kp[p] | 0.1); 
    target += exponential_lpdf(sigma_kp[p] | 0.1); 
  }  
  //target += normal_lpdf(ap|mu_a,2.5);
  target += normal_lpdf(ap|0,10);
}
generated quantities {
  array[P] vector[Tf] kp_p;
  vector[ATfP] mfore; // predicted death rates
  vector[ATP] log_lik;
  array[ATP] real m;
  int pos1 = 1;
  
  for(p in 1:P){
    kp_p[p][1] = cp[p] + kp[p][T] + sigma_kp[p] * normal_rng(0,1);
    for (t in 2:Tf){
      kp_p[p][t] = cp[p] + kp_p[p][t-1] + sigma_kp[p] * normal_rng(0,1);
    }
    for (t in 1:Tf){
      for (i in 1:A) {
      mfore[pos1] = exp(ap[(p-1)*A + i] + bp[p][i] * kp_p[p][t]);
      pos1 += 1;
      }
    } 
    
  }
  for(i in 1:ATP){
    m[i] = exp(logm[i] - loge[i]);
    log_lik[i] = poisson_log_lpmf(death[i]|logm[i]);
  }
  
}
