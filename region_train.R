library(cmdstanr)
#set_cmdstan_path("E:/tools/cmdstan/cmdstan-2.31.0")
library(bayesplot)
library(tidyverse)
library(robustbase)
library(loo)
library(ggplot2)
library(reshape2)

# death <- expos <- rate <- array(0,c(nrow(region_Hokkaido$death$total$train), ncol(region_Hokkaido$death$total$train), 8))
# mu_a <- matrix(0, ncol = 8, nrow = nrow(region_Hokkaido$death$total$train))
# mu_b <- matrix(0, ncol = 8, nrow = ncol(region_Hokkaido$death$total$train))
# 
# for (r in 1:8) {
#   death[,,r] <- get(region[r])$death$total$train
#   expos[,,r] <- get(region[r])$expos$total$train
#   rate[,,r] = get(region[r])$rate$total$train
#   mu_a[,r] <- log(rowMeans(rate[,,r]))
#   mu_b[,r] <- log(colMeans(rate[,,r]))
#   # mu_a[,s] <- log(rowSums(death2[,,s])/rowSums(expos2[,,s]))
#   # mu_b[,s] <- log(colSums(death2[,,s])/colSums(expos2[,,s]))
# }


death <- expos <- rate <- array(0,c(nrow(region_Hokkaido$death$total$train), 23, 8))
test_rate <- test_death <- test_expos <- array(0,c(nrow(region_Hokkaido$death$total$test), ncol(region_Hokkaido$death$total$test), 8))
mu_a <- matrix(0, ncol = 8, nrow = nrow(region_Hokkaido$death$total$train))
mu_b <- matrix(0, ncol = 8, nrow = 23)

for (r in 1:8) {
  death[,,r] <- get(region[r])$death$total$train[,44:66]
  expos[,,r] <- get(region[r])$expos$total$train[,44:66]
  rate[,,r] = get(region[r])$rate$total$train[,44:66]
  test_rate[,,r] = get(region[r])$death$total$test/get(region[r])$expos$total$test
  test_death[,,r] = get(region[r])$death$total$test
  test_expos[,,r] = get(region[r])$expos$total$test
  mu_a[,r] <- log(rowMeans(rate[,,r]))
  mu_b[,r] <- log(colMeans(rate[,,r]))
  # mu_a[,s] <- log(rowSums(death2[,,s])/rowSums(expos2[,,s]))
  # mu_b[,s] <- log(colSums(death2[,,s])/colSums(expos2[,,s]))
}


Tfore = 50

standata <- list(death = as.integer(as.vector(death)), # death
                  expos = as.integer(as.vector(expos)),
                  loge = log(as.integer(as.vector(expos))), # exposure
                  #eq = earthquake,
                  #eq_index = eq_index,
                  A = dim(death)[1],        # number of age categories
                  T = dim(death)[2],           # number of years
                  P = dim(death)[3],
                  AT = dim(death)[1] * dim(death)[2],
                  ATP = dim(death)[1] * dim(death)[2] * dim(death)[3],
                  ATfP = dim(death)[1] * Tfore * dim(death)[3],
                  Tf = Tfore,                              # number of forecast years
                  W = w,
                  n = dim(w)[1],
                  W_n = sum(w) / 2,
                  mu_a = as.vector(mu_a)
)

stanfile <- "stanfile/LC_O.stan"
rdsfile <- "output/region/region_LC_O.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 1000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LC_C.stan"
rdsfile <- "output/region/region_LC_C.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 1000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}
stanfile <- "stanfile/LC_O_2.stan"
rdsfile <- "output/region/region_LC_O_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 1000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}
stanfile <- "stanfile/LC_C.stan"
rdsfile <- "output/region/region_LC_C.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LC_C_2.stan"
rdsfile <- "output/region/region_LC_C_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LL_O.stan"
rdsfile <- "output/region/region_LL_O.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LL_O_2.stan"
rdsfile <- "output/region/region_LL_O_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LL_C.stan"
rdsfile <- "output/region/region_LL_C.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

# stanfile <- "stanfile/LL_C_test.stan"
# rdsfile <- "output/region/region_LL_C_test.RDS"
# if(!file.exists(rdsfile)){
#   model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
#   fit <- model$sample(standata,
#                       chains = 4,
#                       parallel_chains = 4,
#                       iter_warmup = 4000,
#                       iter_sampling = 4000,
#                       threads_per_chain = 6,
#                       # max_treedepth = 15,
#                       # adapt_delta = 0.95,
#                       thin = 4,
#                       refresh = 4000)
#   fit$draws()
#   fit$sampler_diagnostics()
#   saveRDS(fit,rdsfile)
# }

stanfile <- "stanfile/LL_C_2.stan"
rdsfile <- "output/region/region_LL_C_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

#Earthquake
earthquake <- matrix(0, 8, 23)
earthquake[2, 22] = 1
earthquake[5, 6] = 1
eq_index =rep(1, 8)
count = 1
for (i in 1:8) {
  if(max(earthquake[i,] > 0)){
    count = count + 1
    eq_index[i] = count
  }
}

standata <- list(death = as.integer(as.vector(death)), # death
                 expos = as.integer(as.vector(expos)),
                 loge = log(as.integer(as.vector(expos))), # exposure
                 #eq = earthquake,
                 #eq_index = eq_index,
                 A = dim(death)[1],        # number of age categories
                 T = dim(death)[2],           # number of years
                 P = dim(death)[3],
                 AT = dim(death)[1] * dim(death)[2],
                 ATP = dim(death)[1] * dim(death)[2] * dim(death)[3],
                 ATfP = dim(death)[1] * Tfore * dim(death)[3],
                 Tf = Tfore,                              # number of forecast years
                 W = w,
                 n = dim(w)[1],
                 W_n = sum(w) / 2,
                 mu_a = as.vector(mu_a),
                 eq = earthquake,
                 eq_index = eq_index,
                 n_eq = max(eq_index) - 1
)
stanfile <- "stanfile/LC_O_E.stan"
rdsfile <- "output/region/region_LC_O_E.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LC_C_E.stan"
rdsfile <- "output/region/region_LC_C_E.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


stanfile <- "stanfile/LL_O_E.stan"
rdsfile <- "output/region/region_LL_O_E.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LL_C_E.stan"
rdsfile <- "output/region/region_LL_C_E.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

# stanfile <- "stanfile/LL_C_E_test.stan"
# rdsfile <- "output/region/region_LL_C_E_test.RDS"
# if(!file.exists(rdsfile)){
#   model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
#   fit <- model$sample(standata,
#                       chains = 4,
#                       parallel_chains = 4,
#                       iter_warmup = 4000,
#                       iter_sampling = 4000,
#                       threads_per_chain = 6,
#                       # max_treedepth = 15,
#                       # adapt_delta = 0.95,
#                       thin = 4,
#                       refresh = 4000)
#   fit$draws()
#   fit$sampler_diagnostics()
#   saveRDS(fit,rdsfile)
# }
stanfile <- "stanfile/LC_O_E_2.stan"
rdsfile <- "output/region/region_LC_O_E_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LC_C_E_2.stan"
rdsfile <- "output/region/region_LC_C_E_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


stanfile <- "stanfile/LL_O_E_2.stan"
rdsfile <- "output/region/region_LL_O_E_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/LL_C_E_2.stan"
rdsfile <- "output/region/region_LL_C_E_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}





#Earthquake
earthquake <- matrix(0, 8, 23)
earthquake[2, 22] = 1
#earthquake[5, 6] = 1
eq_index =rep(1, 8)
count = 1
for (i in 1:8) {
  if(max(earthquake[i,] > 0)){
    count = count + 1
    eq_index[i] = count
  }
}

standata <- list(death = as.integer(as.vector(death)), # death
                 expos = as.integer(as.vector(expos)),
                 loge = log(as.integer(as.vector(expos))), # exposure
                 #eq = earthquake,
                 #eq_index = eq_index,
                 A = dim(death)[1],        # number of age categories
                 T = dim(death)[2],           # number of years
                 P = dim(death)[3],
                 AT = dim(death)[1] * dim(death)[2],
                 ATP = dim(death)[1] * dim(death)[2] * dim(death)[3],
                 ATfP = dim(death)[1] * Tfore * dim(death)[3],
                 Tf = Tfore,                              # number of forecast years
                 W = w,
                 n = dim(w)[1],
                 W_n = sum(w) / 2,
                 mu_a = as.vector(mu_a),
                 eq = earthquake,
                 eq_index = eq_index,
                 n_eq = max(eq_index) - 1
)
stanfile <- "stanfile/NP/LC_O_E.stan"
rdsfile <- "output/region/A/region_LC_O_EE.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/NP/LC_C_E.stan"
rdsfile <- "output/region/region_LC_C_EE.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


stanfile <- "stanfile/NP/LL_O_E.stan"
rdsfile <- "output/region/region_LL_O_EE.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/NP/LL_C_E.stan"
rdsfile <- "output/region/region_LL_C_EE.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


stanfile <- "stanfile/NP/LC_O_E_2.stan"
rdsfile <- "output/region/region_LC_O_EE_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/NP/LC_C_E_2.stan"
rdsfile <- "output/region/region_LC_C_EE_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


stanfile <- "stanfile/NP/LL_O_E_2.stan"
rdsfile <- "output/region/region_LL_O_EE_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

stanfile <- "stanfile/NP/LL_C_E_2.stan"
rdsfile <- "output/region/region_LL_C_EE_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}


#Earthquake3
earthquake <- matrix(0, 8, 23)
earthquake[2, 22] = 1
earthquake[2, 23] = 1
#earthquake[5, 6] = 1
eq_index =rep(1, 8)
count = 1
for (i in 1:8) {
  if(max(earthquake[i,] > 0)){
    count = count + 1
    eq_index[i] = count
  }
}

standata <- list(death = as.integer(as.vector(death)), # death
                 expos = as.integer(as.vector(expos)),
                 loge = log(as.integer(as.vector(expos))), # exposure
                 #eq = earthquake,
                 #eq_index = eq_index,
                 A = dim(death)[1],        # number of age categories
                 T = dim(death)[2],           # number of years
                 P = dim(death)[3],
                 AT = dim(death)[1] * dim(death)[2],
                 ATP = dim(death)[1] * dim(death)[2] * dim(death)[3],
                 ATfP = dim(death)[1] * Tfore * dim(death)[3],
                 Tf = Tfore,                              # number of forecast years
                 W = w,
                 n = dim(w)[1],
                 W_n = sum(w) / 2,
                 mu_a = as.vector(mu_a),
                 eq = earthquake,
                 eq_index = eq_index,
                 n_eq = max(eq_index) - 1
)

stanfile <- "stanfile/NP/LL_C_E.stan"
rdsfile <- "output/region/region_LL_C_E_test.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata,
                      chains = 4,
                      parallel_chains = 4,
                      iter_warmup = 4000,
                      iter_sampling = 4000,
                      threads_per_chain = 6,
                      # max_treedepth = 15,
                      # adapt_delta = 0.95,
                      thin = 4,
                      refresh = 4000)
  fit$draws()
  fit$sampler_diagnostics()
  saveRDS(fit,rdsfile)
}

# 计算RMSE
rmse <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

# 计算MAE
mae <- function(actual, predicted) {
  mean(abs(predicted - actual))
}
r_squared <- function(actual, predicted) {
  # 计算总平方和（Total Sum of Squares）
  SST <- sum((actual - mean(actual))^2)
  # 计算残差平方和（Residual Sum of Squares）
  SSR <- sum((actual - predicted)^2)
  # 计算R平方
  R2 <- 1 - SSR/SST
  return(R2)
}


getsummarymedian_all <- function(folder = "output/region/"){
  allresult <- list() # 用于存储所有结果的列表
  all_comp <- data.frame() # 初始化空的数据框用于汇总comp
  waic_list <- list() 
  loo_list <- list() 
  result <- list()
  # 获取文件夹中所有.RDS文件
  files <- list.files(path = folder, pattern = "\\.RDS$", full.names = TRUE)
  
  # 遍历所有文件
  for (name in files) {
    # 读取每个文件
    fit <- readRDS(name)
    res <- list() # 存储当前文件结果的列表
    
    # 计算DIC
    log_lik <- fit$draws("log_lik", format = "matrix")
    #D1 <- -2 * sum(log(colMeans(exp(log_lik))))
    #m <- fit$draws("m", format = "matrix")
    #m_hat <- colMedians(m)
    #D2 <- -2 * sum(dpois(death_train, m_hat * expos_train, log = TRUE))
    #DIC <- 2 * D1 - D2
    # # 计算WAIC
    # p_waic <- sum(apply(log_lik, 2, var))
    # WAIC <- D1 + 2 * p_waic
    ### WAIC
    W <- waic(log_lik)
    WAIC <- W$estimates[3,1]
    #result$W <- W
    
    
    # 计算Loo
    LOO <- fit$loo()
    LOOIC <- LOO$estimates[3, 1]
    
    # 未来预测值
    #Japan_mu_fcst <- colMedians(fit$draws("mfore", format = "matrix"))
    
    #result$mf <- Japan_mu_fcst
    
    m <- as.data.frame(fit$summary("m",~quantile(.x, probs = c(0.025, 0.25,0.5,0.75,0.975))))[,2:6]
    mfore <- as.data.frame(fit$summary("mfore",~quantile(.x, probs = c(0.025, 0.25,0.5,0.75,0.975))))[,2:6]
    mf <- list()
    m_i <- list()
    mfore_i <- list()
    Tfore = 50
    for (i in 1:5) {
      m_i[[i]] <- array(m[,i], c(35, 23, 8))
      mfore_i[[i]] <- array(mfore[,i], c(35, Tfore, 8))
      mf[[i]] <- array(0,c(35, 23 + Tfore, 8))
      for (p in 1:8) {
        mf[[i]][,,p] <- cbind(m_i[[i]][,,p], mfore_i[[i]][,,p])
      }
    }
    res$m = m_i
    res$mfore = mfore_i
    res$mf = mf
    
    RMSE = rmse(as.vector(test_rate[,1:9,]), as.vector(mfore_i[[3]][,1:9,]))
    MAE = mae(as.vector(test_rate[,1:9,]), as.vector(mfore_i[[3]][,1:9,]))
    R2 = r_squared(as.vector(test_rate[,1:9,]), as.vector(mfore_i[[3]][,1:9,]))
    # 汇总统计量
    comp <- matrix(c(WAIC, LOOIC, RMSE, MAE, R2), nrow = 1)
    colnames(comp) <- c("WAIC", "LOOIC","RMSE", "MAE", "R2")
    rownames(comp) <- basename(name)
    res$comp <- comp
    
    result[[basename(name)]] <- res
    waic_list[[basename(name)]] <- W
    loo_list[[basename(name)]] <- LOO
    # 将comp添加到all_comp数据框
    all_comp <- rbind(all_comp, comp)
    print(name)
  }

  

  
  
  # 打印汇总的comp表格
  print(all_comp)
  allresult$result <- result
  allresult$waic_list <- waic_list
  allresult$loo_list <- loo_list
  return(allresult)
}

allresult <- getsummarymedian_all("output/region/")
save(allresult, file = "compare/region_A2.RData")

load("compare/region_A.RData")
mfore1 = allresult$result$region_LC_C.RDS$mfore[[3]]
mfore2 = allresult$result$region_LC_C_2.RDS$mfore[[3]]
mfore3 = allresult$result$region_LC_O.RDS$mfore[[3]]
mfore4 = allresult$result$region_LC_O_2.RDS$mfore[[3]]
mfore5 = allresult$result$region_LL_C.RDS$mfore[[3]]
mfore6 = allresult$result$region_LL_C_2.RDS$mfore[[3]]
mfore7 = allresult$result$region_LL_O.RDS$mfore[[3]]
mfore8 = allresult$result$region_LL_O_2.RDS$mfore[[3]]

mfore9 = allresult$result$region_LC_C_E.RDS$mfore[[3]]
mfore10 = allresult$result$region_LC_C_E_2.RDS$mfore[[3]]
mfore11 = allresult$result$region_LC_O_E.RDS$mfore[[3]]
mfore12 = allresult$result$region_LC_O_E_2.RDS$mfore[[3]]
mfore13 = allresult$result$region_LL_C_E.RDS$mfore[[3]]
mfore14 = allresult$result$region_LL_C_E_2.RDS$mfore[[3]]
mfore15 = allresult$result$region_LL_O_E.RDS$mfore[[3]]
mfore16 = allresult$result$region_LL_O_E_2.RDS$mfore[[3]]
# mfore_stacking_NOe = 0.019*mfore3 + 0.491*mfore5 + 0.027*mfore6 + 0.463*mfore7
# mfore_pseudobma_NOe = mfore5
# mfore_stacking_WITHe = 0.14*mfore9 + 0.118*mfore10 + 0.135*mfore11 + 0.122*mfore12 + 0.164*mfore5 + 0.053*mfore6 + 0.007*mfore13 + 0.002*mfore10 + 0.082*mfore7 + 0.043*mfore8 +0.133*mfore16
# mfore_pseudobma_WITHe = 0.006*mfore10 + 0.993*mfore16
# MAE_R_stacking_NOe = MAE_R_pseudobma_NOe = MAE_R_stacking_WITHe = MAE_R_pseudobma_WITHe = 0
# for (i in 1:10) {
#   MAE_R_stacking_NOe[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore_stacking_NOe[,i,]))
#   MAE_R_pseudobma_NOe[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore_pseudobma_NOe[,i,]))
#   MAE_R_stacking_WITHe[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore_stacking_WITHe[,i,]))
#   MAE_R_pseudobma_WITHe[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore_pseudobma_WITHe[,i,]))
#   
# }
# MAE_R_stacking_NOe
# MAE_R_pseudobma_NOe
# MAE_R_stacking_WITHe
# MAE_R_pseudobma_WITHe
MAE1 = 0
MAE2 = 0
for (i in 1:10) {
  MAE1[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore5[,i,]))
  MAE2[i] = mae(as.vector(test_rate[,i,]), as.vector(mfore16[,i,]))
#   cumulate_MAE[i] = mae(as.vector(test_rate[,1:i,]), as.vector(mfore5[,1:i,]))
 }
# waics <- 0
# for (i in 1:length(allresult$waic_list)) {
#   waics[i] = allresult$waic_list[[i]]$estimates["elpd_waic", 1]
# }
# max_waic <- max(waics)
# adjusted_waics <- waics - max_waic
# 对调整后的值进行指数运算并计算权重
# exp_adjusted_waics <- exp(adjusted_waics)
# waic_wts <- exp_adjusted_waics / sum(exp_adjusted_waics)
# lpd_point <- allresult$loo_list$region_LC_C.RDS$pointwise[,"elpd_loo"]
# 
# for (i in 2:length(allresult$loo_list)) {
#   lpd_point = cbind(lpd_point, allresult$loo_list[[i]]$pointwise[,"elpd_loo"])
# }
# pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
# pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
# stacking_wts <- stacking_weights(lpd_point)
# round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2)


library(dplyr)

fit <- readRDS(files[1])

mfore <- as.array(fit$draws("mfore"))

dval <- as.integer(as.vector(test_death))
eval <- as.integer(as.vector(test_expos))

calc_log_lik <- function(mfore, death, expos) {
  n_obs <- length(death)
  n_iter <- dim(mfore)[1]
  n_chain <- dim(mfore)[2]
  
  log_lik <- array(NA, dim = c(n_iter, n_chain, n_obs))
  
  for (i in 1:n_obs) {
    log_lik[, , i] <- dpois(death[i], lambda = mfore[, , i]*expos[i], log = TRUE)
  }
  
  return(log_lik)
}


#log_lik2 <- calc_log_lik(mfore, dval, eval)

# calc_weights <- function(path = "output/region/A/", dval, eval) {
#   log_sum_exp <- function(u) {
#     max_u <- max(u)
#     a <- 0
#     for (n in 1:length(u)) {
#       a <- a + exp(u[n] - max_u)
#     }
#     return(max_u + log(a))
#   }
#   
#   files <- list.files(path = path, pattern = "\\.RDS$", full.names = TRUE)
#   
#   # 创建进度条
#   pb <- txtProgressBar(min = 0, max = length(files), style = 3)
#   
#   log_lik_list <- lapply(seq_along(files), function(i) {
#     file <- files[i]
#     fit <- readRDS(file)
#     mfore <- as.array(fit$draws("mfore"))
#     
#     # 更新进度条
#     setTxtProgressBar(pb, i)
#     
#     calc_log_lik(mfore, dval, eval)
#   })
#   
#   # 关闭进度条
#   close(pb)
#   
#   lpd <- lapply(log_lik_list, function(log_lik2) {
#     apply(log_lik2, 3, log_sum_exp) - log(dim(log_lik2)[1])
#   })
#   
#   lpd_point <- simplify2array(lpd)
#   stacking <- loo::stacking_weights(lpd_point)
#   pseudobma <- loo::pseudobma_weights(lpd_point, BB = FALSE)
#   
#   return(cbind(files, round(cbind(stacking, pseudobma), 3)))
# }
# # 计算模型权重
# weights <- calc_weights("output/region/A", dval, eval)
# weights
