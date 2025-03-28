functions {
  real partial_sum_lpdf(array[] real logm,
                        int start, int end,
                        array[] int death) {
    return poisson_log_lupmf(death[start:end] | logm);
  }
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Scale parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    array[,] int W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (-2 * n * log(tau)
                    + sum(ldet_terms)
                    - (1/tau^2) * (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data {
  int<lower=0> T;                  // number of years
  int<lower=0> A;                  // number of age categories
  int<lower=0> P;                  // number of age categories
  int<lower=0> Tf;                 // number of forecast years
  int<lower=0> AT;
  int<lower=0> ATP;
  int<lower=0> ATfP;
  array[ATP] real loge;        // log exposures
  array[ATP] int<lower=0> death;   // deaths
  int<lower=0> n; 
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
}
transformed data {
  int grainsize = 1;
  array[W_n, 2] int W_sparse;   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n) D_sparse[i] = sum(W[i]);
  {
    vector[n] invsqrtD;  
    for (i in 1:n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters {
  simplex[A] b; 
  real c;
  vector[T-1] ks;                 // vector of k_t differences
  real<lower = 0> sigma_k;          // standard deviation of the random walk
  
  array[P] simplex[A] bp; 
  array[P] real cp;
  array[P] vector[T-1] kps;                 // vector of k_t differences
  array[P] real<lower = 0> sigma_kp;          // standard deviation of the random walk
  
  array[A] vector[n] theta;
  vector<lower = 0>[A] tau;
  vector<lower = 0, upper = 1>[A] alpha;
}
transformed parameters  {
  vector[T] k;
  array[P] vector[T] kp; 
  array[ATP] real logm;
  
  k[1] = 0;
  //k[2] = -sum(ks[1:(T-2)]);
  k[2:T] = ks;
  for(p in 1:P){
    kp[p][1] = 0;
    //kp[p][2] = -sum(kps[p][1:(T-2)]);
    kp[p][2:T] = kps[p];
  }
  
  for(p in 1:P){
    for (j in 1:T){
      for (i in 1:A) {
        logm[AT*(p-1) + A * (j-1) + i] = loge[AT*(p-1) + A * (j-1) + i] + 
                                         b[i] * k[j] + 
                                         bp[p][i] * kp[p][j] +
                                         theta[i][p];
      }
    }
  }
  
  
}
model {
  target += reduce_sum(partial_sum_lpdf, logm, grainsize, death);
  //target += poisson_log_lpmf(death| logm);   
  target += normal_lpdf(ks[1]| c, sigma_k);
  target += normal_lpdf(ks[2:(T-1)]| c + ks[1:(T-2)], sigma_k);
  target += dirichlet_lpdf(b|rep_vector(1, A));
  target += normal_lpdf(c|0, sqrt(10));
  //target += student_t_lupdf(sigma_k|4, 0, 1);
  target += exponential_lpdf(sigma_k | 0.1); 
  for(p in 1:P){
    target += normal_lpdf(kps[p][1]| cp[p], sigma_kp[p]);
    target += normal_lpdf(kps[p][2:(T-1)]| cp[p] + kps[p][1:(T-2)], sigma_kp[p]);
    target += dirichlet_lpdf(bp[p]|rep_vector(1, A));
    target += normal_lpdf(cp[p]|0, sqrt(10));
    //target += normal_lpdf(sigma_k|2, 1);
    //target += student_t_lupdf(sigma_kp[p]|4, 0, 1);
    target += exponential_lpdf(sigma_kp[p] | 0.1); 
  }  
  for(i in 1:A){
    target += sparse_car_lpdf(theta[i] |tau[i], alpha[i], W_sparse, D_sparse, lambda, n, W_n);
    target += std_normal_lpdf(tau[i]);
  }
  
}
generated quantities {
  vector[Tf] k_p;
  array[P] vector[Tf] kp_p;
  vector[ATfP] mfore; // predicted death rates
  vector[ATP] log_lik;
  array[ATP] real m;
  int pos1 = 1;
  
  k_p[1] = c + k[T] + sigma_k * normal_rng(0,1);
  for (t in 2:Tf){
    k_p[t] = c + k_p[t-1] + sigma_k * normal_rng(0,1);
  }
  for(p in 1:P){
    kp_p[p][1] = cp[p] + kp[p][T] + sigma_kp[p] * normal_rng(0,1);
    for (t in 2:Tf){
      kp_p[p][t] = cp[p] + kp_p[p][t-1] + sigma_kp[p] * normal_rng(0,1);
    }
    for (t in 1:Tf){
      for (i in 1:A) {
      mfore[pos1] = exp(b[i] * k_p[t] +
                        bp[p][i] * kp_p[p][t] +
                        theta[i][p]);
      pos1 += 1;
      }
    } 
    
  }
  for(i in 1:ATP){
    m[i] = exp(logm[i] - loge[i]);
    log_lik[i] = poisson_log_lpmf(death[i]|logm[i]);
  }
  
}
