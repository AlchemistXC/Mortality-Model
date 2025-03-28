load("compare/region_A2.RData")
allresultR = allresult
load("compare/prefecture_A2.RData")
allresultP = allresult
MAE_LLC_t_R = MAE_LLC_e_R = MAE_LCCE_t_R = MAE_LCCE_e_R = MAE_LLC_t_P = MAE_LLC_e_P = MAE_LCCE_t_P = MAE_LCCE_e_P = rep(0, 10)
for (i in 1:10) {
  MAE_LLC_t_R[i] = mae(as.vector(test_rate[,i,]), as.vector(allresultR$result$region_LL_C.RDS$mfore[[3]][,i,]))
  MAE_LLC_e_R[i] = mae(as.vector(test_rate[,i,]), as.vector(allresultR$result$region_LL_C_2.RDS$mfore[[3]][,i,]))
  MAE_LCCE_t_R[i] = mae(as.vector(test_rate[,i,]), as.vector(allresultR$result$region_LC_C_E.RDS$mfore[[3]][,i,]))
  MAE_LCCE_e_R[i] = mae(as.vector(test_rate[,i,]), as.vector(allresultR$result$region_LC_C_E_2.RDS$mfore[[3]][,i,]))
  MAE_LLC_t_P[i] = mae(as.vector(test_rate2[,i,]), as.vector(allresultP$result$prefecture_LL_C.RDS$mfore[[3]][,i,]))
  MAE_LLC_e_P[i] = mae(as.vector(test_rate2[,i,]), as.vector(allresultP$result$prefecture_LL_C_2.RDS$mfore[[3]][,i,]))
  MAE_LCCE_t_P[i] = mae(as.vector(test_rate2[,i,]), as.vector(allresultP$result$prefecture_LC_C_E.RDS$mfore[[3]][,i,]))
  MAE_LCCE_e_P[i] = mae(as.vector(test_rate2[,i,]), as.vector(allresultP$result$prefecture_LC_C_E_2.RDS$mfore[[3]][,i,]))
}


# 创建年份向量
years <- 2013:2022

# 创建数据框
results_df <- data.frame(
  Year = years,
  MAE_LLC_t_R = MAE_LLC_t_R,
  MAE_LLC_e_R = MAE_LLC_e_R,
  MAE_LCCE_t_R = MAE_LCCE_t_R,
  MAE_LCCE_e_R = MAE_LCCE_e_R,
  MAE_LLC_t_P = MAE_LLC_t_P,
  MAE_LLC_e_P = MAE_LLC_e_P,
  MAE_LCCE_t_P = MAE_LCCE_t_P,
  MAE_LCCE_e_P = MAE_LCCE_e_P
)

# 写入CSV文件
write.csv(results_df, file = "results.csv", row.names = TRUE)
