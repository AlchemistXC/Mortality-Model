# 加载所需包
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(reshape2)

raw_data <- Japan$rate$total$raw[,44:76]
write.csv(log(raw_data), "japan_rate.csv")
long_data = melt(raw_data, varnames = c("Age", "Year"), value.name = "rate")
long_data$log_rate = log(long_data$rate)
write.csv(t(log(raw_data)), "japan_rate2.csv")
ggplot(long_data, aes(x = Year, y = Age, fill = log_rate)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(long_data$log_rate)) +
  labs(x = "Year", y = "Age") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(), # 移除网格线
        panel.background = element_blank(), # 移除背景
        legend.position = "right")
ggsave(file="plot/Log_Mortality_rate_of_Japan.pdf", width=15, height=15)
ggsave(file="plot/Log_Mortality_rate_of_Japan.jpg", width=15, height=15)


raw_data
ggplot(long_data, aes(x=Year, y=log_rate, group=Age, color=as.numeric(Age))) +
  geom_line(linewidth = 0.7) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray80", linewidth = 0.2),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.1),
        legend.position = "none",
        panel.ontop = TRUE) +  # 设置panel.ontop为TRUE
  labs(x="Year", y="Logarithmic mortality rate") +
  scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "red")) +
  scale_x_continuous(breaks = seq(1990, 2020, 5), expand = c(0, 0)) +  # 设置x轴的expand为(0, 0)
  scale_y_continuous(breaks = seq(-10, 2.5, 0.5), expand = c(0, 0))  # 设置y轴的expand为(0, 0)
ggsave(file="plot/Line_Log_Mortality_rate_of_Japan.pdf", width=15, height=15)
ggsave(file="plot/Line_Log_Mortality_rate_of_Japan.jpg", width=15, height=15)
