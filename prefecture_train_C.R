library(cmdstanr)
#set_cmdstan_path("E:/tools/cmdstan/cmdstan-2.31.0")
library(bayesplot)
library(tidyverse)
library(robustbase)
library(loo)
library(ggplot2)
library(reshape2)
# death2 <- expos2 <- array(0,c(nrow(Hokkaido$death$total$train), ncol(Hokkaido$death$total$train), 46))
# rate2 <- array(0, c(nrow(Hokkaido$death$total$train), ncol(Hokkaido$death$total$train), 46))
# real_rate2 <- array(0, c(nrow(Hokkaido$death$total$train), ncol(Hokkaido$death$total$train), 46))
# mu_a <- matrix(0, ncol = 46, nrow = nrow(Hokkaido$death$total$train))
# mu_b <- matrix(0, ncol = 46, nrow = ncol(Hokkaido$death$total$train))
# for (s in 1:46) {
#   death2[,,s] <- get(state[s])$death$total$train
#   expos2[,,s] <- get(state[s])$expos$total$train
#   # rate2[,,s] = r = get(state[s])$rate$total$raw[2:51,]
#   real_rate2[,,s] = r = death2[,,s]/expos2[,,s]
#   mu_a[,s] <- log(rowMeans(r))
#   mu_b[,s] <- log(colMeans(r))
#   # mu_a[,s] <- log(rowSums(death2[,,s])/rowSums(expos2[,,s]))
#   # mu_b[,s] <- log(colSums(death2[,,s])/colSums(expos2[,,s]))
# }

death2 <- expos2 <- array(0,c(nrow(Hokkaido$death$total$raw), 29, 46))
rate2 <- array(0, c(nrow(Hokkaido$death$total$raw), 29, 46))
real_rate2 <- array(0, c(nrow(Hokkaido$death$total$raw), 29, 46))
test_rate2 <- array(0, c(nrow(Hokkaido$death$total$raw), 10, 46))
mu_a <- matrix(0, ncol = 46, nrow = nrow(Hokkaido$death$total$raw))
mu_b <- matrix(0, ncol = 46, nrow = 29)
for (s in 1:46) {
  death2[,,s] <- get(state[s])$death$total$raw[,44:72]
  expos2[,,s] <- get(state[s])$expos$total$raw[,44:72]
  # rate2[,,s] = r = get(state[s])$rate$total$raw[2:51,]
  real_rate2[,,s] = r = death2[,,s]/expos2[,,s]
  test_rate2[,,s] = get(state[s])$death$total$test/get(state[s])$expos$total$test
  mu_a[,s] <- log(rowMeans(r))
  mu_b[,s] <- log(colMeans(r))
  # mu_a[,s] <- log(rowSums(death2[,,s])/rowSums(expos2[,,s]))
  # mu_b[,s] <- log(colSums(death2[,,s])/colSums(expos2[,,s]))
}

Tfore = 10

standata2 <- list(death = as.integer(as.vector(death2)), # death
                  expos = as.integer(as.vector(expos2)),
                  loge = log(as.integer(as.vector(expos2))), # exposure
                  #eq = earthquake,
                  #eq_index = eq_index,
                  A = dim(death2)[1],        # number of age categories
                  T = dim(death2)[2],           # number of years
                  P = dim(death2)[3],
                  AT = dim(death2)[1] * dim(death2)[2],
                  ATP = dim(death2)[1] * dim(death2)[2] * dim(death2)[3],
                  ATfP = dim(death2)[1] * Tfore * dim(death2)[3],
                  Tf = Tfore,                              # number of forecast years
                  W = w,
                  n = dim(w)[1],
                  W_n = sum(w) / 2,
                  mu_a = as.vector(mu_a)
)

stanfile <- "stanfile/LC_O.stan"
rdsfile <- "output/prefecture/prefecture_LC_O.RDS"

if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
rdsfile <- "output/prefecture/prefecture_LC_C.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
rdsfile <- "output/prefecture/prefecture_LC_O_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
# stanfile <- "stanfile/LC_C.stan"
# rdsfile <- "output/prefecture/prefecture_LC_C.RDS"
# if(!file.exists(rdsfile)){
#   model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
#   fit <- model$sample(standata2,
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

stanfile <- "stanfile/LC_C_2.stan"
rdsfile <- "output/prefecture/prefecture_LC_C_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
rdsfile <- "output/prefecture/prefecture_LL_O.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
rdsfile <- "output/prefecture/prefecture_LL_O_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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
rdsfile <- "output/prefecture/prefecture_LL_C.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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

stanfile <- "stanfile/LL_C_2.stan"
rdsfile <- "output/prefecture/prefecture_LL_C_2.RDS"
if(!file.exists(rdsfile)){
  model = cmdstan_model(stanfile, cpp_options = list(stan_threads = TRUE))
  fit <- model$sample(standata2,
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


# eq
earthquake <- matrix(0, 46, 29)
earthquake[3, 22] = 1
earthquake[4, 22] = 1
earthquake[7, 22] = 1
earthquake[28, 6] = 1
eq_index =rep(1, 46)
count = 1
for (i in 1:8) {
  if(max(earthquake[i,] > 0)){
    count = count + 1
    eq_index[i] = count
  }
}

standata <- list(death = as.integer(as.vector(death2)), # death
                 expos = as.integer(as.vector(expos2)),
                 loge = log(as.integer(as.vector(expos2))), # exposure
                 #eq = earthquake,
                 #eq_index = eq_index,
                 A = dim(death2)[1],        # number of age categories
                 T = dim(death2)[2],           # number of years
                 P = dim(death2)[3],
                 AT = dim(death2)[1] * dim(death2)[2],
                 ATP = dim(death2)[1] * dim(death2)[2] * dim(death2)[3],
                 ATfP = dim(death2)[1] * Tfore * dim(death2)[3],
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
rdsfile <- "output/prefecture/prefecture_LC_O_E.RDS"
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
rdsfile <- "output/prefecture/prefecture_LC_C_E.RDS"
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
rdsfile <- "output/prefecture/prefecture_LL_O_E.RDS"
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
rdsfile <- "output/prefecture/prefecture_LL_C_E.RDS"
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


stanfile <- "stanfile/LC_O_E_2.stan"
rdsfile <- "output/prefecture/C/prefecture_LC_O_E_2.RDS"
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
rdsfile <- "output/prefecture/C/prefecture_LC_C_E_2.RDS"
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
rdsfile <- "output/prefecture/C/prefecture_LL_O_E_2.RDS"
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
rdsfile <- "output/prefecture/C/prefecture_LL_C_E_2.RDS"
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

getsummarymedian_all <- function(folder = "output/prefecture/C"){
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
    res <- list()
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
    #res$W <- W
    
    
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
    Tfore = 10
    for (i in 1:5) {
      m_i[[i]] <- array(m[,i], c(35, 29, 46))
      mfore_i[[i]] <- array(mfore[,i], c(35, Tfore, 46))
      mf[[i]] <- array(0,c(35, 29 + Tfore, 46))
      for (p in 1:46) {
        mf[[i]][,,p] <- cbind(m_i[[i]][,,p], mfore_i[[i]][,,p])
      }
    }
    res$m = m_i
    res$mfore = mfore_i
    res$mf = mf
    # 汇总统计量
    comp <- matrix(c(WAIC, LOOIC), nrow = 1)
    colnames(comp) <- c("WAIC", "LOOIC")
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

allresult <- getsummarymedian_all("output/prefecture/C")
save(allresult, file = "compare/prefecture_C.RData")


