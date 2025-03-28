rate = 0
for (i in 1:10) {
  rate[i] = sum(Japan$death$total$test[,i])/sum(Japan$expos$total$test[,i])
}
rate

Japan$rate$total$test


# 创建数据框
# df <- data.frame(
#   Year = 2013:2022,
#   MAE_P_stacking_NOe = c(0.001720412, 0.001331806, 0.001637874, 0.001792806, 0.002350314, 0.002155284, 0.002155766, 0.001252205, 0.002098121, 0.005352497),
#   MAE_P_pseudobma_NOe = c(0.001722878, 0.001332322, 0.001640014, 0.001802982, 0.002368736, 0.002172120, 0.002172874, 0.001255269, 0.002113260, 0.005373546),
#   MAE_P_stacking_WITHe = c(0.001731427, 0.001330565, 0.001653617, 0.001796825, 0.002363960, 0.002162152, 0.002155496, 0.001263184, 0.002071811, 0.005314158),
#   MAE_P_pseudobma_WITHe = c(0.001722878, 0.001332322, 0.001640014, 0.001802982, 0.002368736, 0.002172120, 0.002172874, 0.001255269, 0.002113260, 0.005373546),
#   MAE_R_stacking_NOe = c(0.0012394096, 0.0009250172, 0.0013609407, 0.0014412123, 0.0019539383, 0.0018428127, 0.0018905624, 0.0007538173, 0.0017976679, 0.0050332659),
#   MAE_R_pseudobma_NOe = c(0.0011988561, 0.0009118584, 0.0013191369, 0.0014193602, 0.0019407080, 0.0018182646, 0.0018777221, 0.0007686228, 0.0018026063, 0.0050400892),
#   MAE_R_stacking_WITHe = c(0.0013578442, 0.0009758200, 0.0014626187, 0.0015691981, 0.0020588290, 0.0019418445, 0.0019795551, 0.0007800166, 0.0018646459, 0.0050879747),
#   MAE_R_pseudobma_WITHe = c(0.0012394096, 0.0009250172, 0.0013609407, 0.0014412123, 0.0019539383, 0.0018428127, 0.0018905624, 0.0007538173, 0.0017976679, 0.0050332659)
# )
df <- data.frame(
  Year = 2013:2022,
  llc_r = MAE1,
  lloe_r = MAE2,
  llc_p = MAE3,
  lcce_p = MAE4
)
write.csv(df, "csv/covid19_mae.csv")
library(ggplot2)
library(tidyr)
library(dplyr)

# df_long <- df %>%
#   pivot_longer(cols = -Year, names_to = "Model", values_to = "MAE") %>%
#   mutate(Level = ifelse(grepl("_P_", Model), "Prefecture", "Region"),
#          Method = case_when(
#            grepl("pseudobma_NOe", Model) ~ "PBMA_NoE",
#            grepl("stacking_NOe", Model) ~ "STK_NoE",
#            grepl("stacking_WITHe", Model) ~ "STK_WE",
#            grepl("pseudobma_WITHe", Model) ~ "PBMA_WE"
#          ))

df_long <- df %>%
  pivot_longer(cols = -Year, names_to = "Model", values_to = "MAE")
ggplot(df_long, aes(x = Year, y = MAE, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Prefecture" = "blue", "Region" = "red")) +
  labs(x = "Year", y = "MAE", color = "Level", linetype = "Method") +
  scale_x_continuous(breaks = seq(2013, 2023, 2), expand = c(0, 0)) + 
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file="plot/Mae_COVID-19.pdf", width=15, height=10)
ggsave(file="plot/Mae_COVID-19.jpg", width=15, height=10)
