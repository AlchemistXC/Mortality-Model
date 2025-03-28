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


death <- expos <- rate <- array(0,c(nrow(region_Hokkaido$death$total$raw), 29, 8))
test_rate <- array(0,c(nrow(region_Hokkaido$death$total$test), ncol(region_Hokkaido$death$total$test), 8))
mu_a <- matrix(0, ncol = 8, nrow = nrow(region_Hokkaido$death$total$train))
mu_b <- matrix(0, ncol = 8, nrow = 29)

for (r in 1:8) {
  death[,,r] <- get(region[r])$death$total$raw[,44:72]
  expos[,,r] <- get(region[r])$expos$total$raw[,44:72]
  rate[,,r] = get(region[r])$rate$total$raw[,44:72]
  test_rate[,,r] = get(region[r])$death$total$test/get(region[r])$expos$total$test
  mu_a[,r] <- log(rowMeans(rate[,,r]))
  mu_b[,r] <- log(colMeans(rate[,,r]))
  # mu_a[,s] <- log(rowSums(death2[,,s])/rowSums(expos2[,,s]))
  # mu_b[,s] <- log(colSums(death2[,,s])/colSums(expos2[,,s]))
}


Tfore = 10

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
earthquake <- matrix(0, 8, 29)
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



stanfile <- "stanfile/LC_O_E_2.stan"
rdsfile <- "output/region/C/region_LC_O_E_2.RDS"
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
rdsfile <- "output/region/C/region_LC_C_E_2.RDS"
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
rdsfile <- "output/region/C/region_LL_O_E_2.RDS"
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
rdsfile <- "output/region/C/region_LL_C_E_2.RDS"
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
    Tfore = 10
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
    
    # 汇总统计量
    comp <- matrix(c(WAIC, LOOIC), nrow = 1)
    colnames(comp) <- c("WAIC", "LOOIC")
    rownames(comp) <- basename(name)
    res$comp <- comp
    
    result[[basename(name)]] <- res
    # 将comp添加到all_comp数据框
    all_comp <- rbind(all_comp, comp)
    print(name)
  }
  
  
  
  
  
  print(all_comp)
  allresult$result <- result
  allresult$waic_list <- waic_list
  allresult$loo_list <- loo_list
  return(allresult)
}

allresult <- getsummarymedian_all("output/region/C")
save(allresult, file = "compare/region_C.RData")
