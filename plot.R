library(ggplot2)
library(reshape2)
library(viridis)

allrate1 = matrix(0, ncol = 33, nrow = 47)
allrate2 = matrix(0, ncol = 33, nrow = 46)
allrate3 = matrix(0, ncol = 33, nrow = 46)
allrate4 = matrix(0, ncol = 33, nrow = 46)
allrate = matrix(0, ncol = 33, nrow = 46)
for (i in 2:47) {
  deaths1 = colSums(get(state[i-1])$death$total$raw[,44:76])
  exp1 = colSums(get(state[i-1])$expos$total$raw[,44:76])
  allrate1[i,] = deaths1 / exp1
  # 
  # deaths2 = colSums(get(state[i])$death$total$raw[14:23,44:76])
  # exp2 = colSums(get(state[i])$expos$total$raw[14:23,44:76])
  # allrate2[i,] = deaths2 / exp2
  # 
  # deaths3 = colSums(get(state[i])$death$total$raw[24:35,44:76])
  # exp3 = colSums(get(state[i])$expos$total$raw[24:35,44:76])
  # allrate3[i,] = deaths3 / exp3
  # 
  # deaths4 = get(state[i])$death$total$raw[35,44:76]
  # exp4 = get(state[i])$expos$total$raw[35,44:76]
  # allrate4[i,] = deaths4 / exp4
  # 
  # deaths = colSums(get(state[i])$death$total$raw[,44:76])
  # exp = colSums(get(state[i])$expos$total$raw[,44:76])
  # allrate[i,] = deaths / exp
  # #allrate[i,] = colSums(get(state[i])$rate$total$raw)/46
}
allrate1[1,] = colSums(Japan$death$total$raw[,44:76])/colSums(Japan$expos$total$raw[,44:76])
colnames(allrate1) = colnames(allrate2) = colnames(allrate3) = colnames(allrate4) = colnames(allrate) = 1990:2022  # 假设您的列是从2000年到2020年
rownames(allrate1) = c("Japan",state)
# rownames(allrate2) = rownames(allrate3) = rownames(allrate4) = rownames(allrate) = state # 'state' 应该是包含所有省份名的向量
# rownames(allrate1) = state[1:47]

log(t(allrate1))
dim(allrate1)
write.csv(log(t(allrate1)), "Japan_rate.csv")
long_data1 = melt(allrate1, varnames = c("Prefecture", "Year"), value.name = "rate")
long_data2 = melt(allrate2, varnames = c("Prefecture", "Year"), value.name = "rate")
long_data3 = melt(allrate3, varnames = c("Prefecture", "Year"), value.name = "rate")
long_data4 = melt(allrate4, varnames = c("Prefecture", "Year"), value.name = "rate")
long_data = melt(allrate, varnames = c("Prefecture", "Year"), value.name = "rate")

# ggplot(long_data, aes(x=年份, y=省份, fill=死亡率)) + 
#   geom_tile() + 
#   scale_fill_gradient(low="blue", high="red") +
#   theme_minimal() +
#   labs(title="省份按年份的死亡率", x="年份", y="省份", fill="死亡率")

long_data1$Log_rate = log(long_data1$rate)
long_data2$Log_rate = log(long_data2$rate)
long_data3$Log_rate = log(long_data3$rate)
long_data4$Log_rate = log(long_data4$rate)
long_data$Log_rate = log(long_data$rate)
#----

# 为离散变量着色  
ggplot(long_data1, aes(x=Year, y=Log_rate, group=Prefecture)) +
  geom_line(data = subset(long_data1, Prefecture != "Japan"), 
            aes(color=Prefecture), size = 0.8) +
  geom_line(data = subset(long_data1, Prefecture == "Japan"),
            color = "red", size = 1.5) +
  scale_color_viridis(discrete = TRUE, option = "D", na.value = "red") +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(1990, 2020, 5), expand = c(0, 0)) +  # 设置x轴的expand为(0, 0)
  scale_y_continuous(breaks = seq(-10, 2.5, 0.5), expand = c(0, 0)) +
  labs(x = "Year", y = "Logarithmic mortality rate", 
       color = "Prefecture")


ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve.pdf", width=15, height=10)
ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve.jpg", width=15, height=10)

#----
ggplot(long_data2, aes(x=Year, y=Log_rate, group=Prefecture, color=Prefecture)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Prefectural mortality rate logarithmic curve (Age:34-66)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:34-66).pdf", width=15, height=15)
ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:34-66).jpg", width=15, height=15)

#----
ggplot(long_data3, aes(x=Year, y=Log_rate, group=Prefecture, color=Prefecture)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Prefectural mortality rate logarithmic curve (Age:67-100+)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:67-100+).pdf", width=15, height=15)
ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:67-100+).jpg", width=15, height=15)
#----
ggplot(long_data4, aes(x=Year, y=Log_rate, group=Prefecture, color=Prefecture)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Prefectural mortality rate logarithmic curve (Age:67-100+)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:67-100+).pdf", width=15, height=15)
ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:67-100+).jpg", width=15, height=15)

#----
ggplot(long_data, aes(x=Year, y=Log_rate, group=Prefecture, color=Prefecture)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Prefectural mortality rate logarithmic curve (Age:0-100+)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:0-100+).pdf", width=15, height=15)
ggsave(file="plot/Prefectural_mortality_rate_logarithmic_curve_(Age:0-100+).jpg", width=15, height=15)


#---------------------------------------------------------------------------------------------------------------------------------------

allrate1 = matrix(0, ncol = 33, nrow = 9)
allrate2 = matrix(0, ncol = 33, nrow = 8)
allrate3 = matrix(0, ncol = 33, nrow = 8)
allrate = matrix(0, ncol = 33, nrow = 8)
for (i in 1:8) {
  deaths1 = colSums(get(region[i])$death$total$raw[,44:76])
  exp1 = colSums(get(region[i])$expos$total$raw[,44:76])
  allrate1[i+1,] = deaths1 / exp1
  
  deaths2 = colSums(get(region[i])$death$total$raw[14:23,44:76])
  exp2 = colSums(get(region[i])$expos$total$raw[14:23,44:76])
  allrate2[i,] = deaths2 / exp2
  
  deaths3 = colSums(get(region[i])$death$total$raw[24:35,44:76])
  exp3 = colSums(get(region[i])$expos$total$raw[24:35,44:76])
  allrate3[i,] = deaths3 / exp3
  
  deaths = colSums(get(region[i])$death$total$raw[,44:76])
  exp = colSums(get(region[i])$expos$total$raw[,44:76])
  allrate[i,] = deaths / exp
}
allrate1[1,] = colSums(Japan$death$total$raw[,44:76])/colSums(Japan$expos$total$raw[,44:76])
colnames(allrate1) = colnames(allrate2) = colnames(allrate3) = colnames(allrate) = 1990:2022  # 假设您的列是从2000年到2020年
rownames(allrate2) = rownames(allrate3) = rownames(allrate) = region[1:8] # 'state' 应该是包含所有省份名的向量
rownames(allrate1) = c("Japan", "region_Hokkaido", "region_Tohoku", "region_Kanto", "region_Chubu", "region_Kinki", "region_Chugoku", "region_Shikoku", "region_Kyushu")
t(allrate1)
write.csv(log(t(allrate1)), "Japan_region_rate.csv")
long_data1 = melt(allrate1, varnames = c("Region", "Year"), value.name = "rate")
long_data2 = melt(allrate2, varnames = c("Region", "Year"), value.name = "rate")
long_data3 = melt(allrate3, varnames = c("Region", "Year"), value.name = "rate")
long_data = melt(allrate, varnames = c("Region", "Year"), value.name = "rate")

# ggplot(long_data, aes(x=年份, y=省份, fill=死亡率)) + 
#   geom_tile() + 
#   scale_fill_gradient(low="blue", high="red") +
#   theme_minimal() +
#   labs(title="省份按年份的死亡率", x="年份", y="省份", fill="死亡率")

long_data1$Log_rate = log(long_data1$rate)
long_data2$Log_rate = log(long_data2$rate)
long_data3$Log_rate = log(long_data3$rate)
long_data$Log_rate = log(long_data$rate)
#----
ggplot(long_data1, aes(x=Year, y=Log_rate, group=Region)) +
  geom_line(data = subset(long_data1, Region != "Japan"), 
            aes(color=Region), size = 0.8) +
  geom_line(data = subset(long_data1, Region == "Japan"),
            color = "red", size = 1.5) +
  scale_color_viridis(discrete = TRUE, option = "D", na.value = "red") +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(1990, 2020, 5), expand = c(0, 0)) +  # 设置x轴的expand为(0, 0)
  scale_y_continuous(breaks = seq(-10, 2.5, 0.5), expand = c(0, 0)) +
  labs(x = "Year", y = "Logarithmic mortality rate", 
       color = "Prefecture")

ggsave(file="plot/Regional_mortality_rate_logarithmic_curve.pdf", width=15, height=10)
ggsave(file="plot/Regional_mortality_rate_logarithmic_curve.jpg", width=15, height=10)

#----
ggplot(long_data2, aes(x=Year, y=Log_rate, group=Region, color=Region)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Regional mortality rate logarithmic curve (Age:34-66)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_(Age:34-66).pdf", width=15, height=15)
ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_(Age:34-66).jpg", width=15, height=15)

#----
ggplot(long_data3, aes(x=Year, y=Log_rate, group=Region, color=Region)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Regional mortality rate logarithmic curve (Age:67-100+)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_(Age:67-100+).pdf", width=15, height=15)
ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_(Age:67-100+).jpg", width=15, height=15)

#----
ggplot(long_data, aes(x=Year, y=Log_rate, group=Region, color=Region)) + 
  geom_line() +
  theme_minimal() +
  labs(title="Regional mortality rate logarithmic curve (Age:0-100+)", x="Year", y="Logarithmic mortality rate")

ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_Age:0-100+.pdf", width=15, height=15)
ggsave(file="plot/Regional_mortality_rate_logarithmic_curve_Age:0-100+.jpg", width=15, height=15)

colSums(Japan$expos$total$raw[22:35,44:76])/colSums(Japan$expos$total$raw[,44:76])
