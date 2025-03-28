library(ggplot2)
library(tidyr)

death_rate_df <- as.data.frame(log(Fukushima$rate$total$raw))
death_rate_df <- as.data.frame(Hyogo$rate$total$raw)
# 给数据框添加年龄作为新的列
death_rate_df$age <- rownames(Hyogo$rate$total$raw)

# 使用gather()之前，需要确保数据是数据框格式
# 将数据框从宽格式转换为长格式
death_rate_long <- gather(death_rate_df, year, death_rate, -age, factor_key=TRUE)

# 绘制图像
ggplot(death_rate_long, aes(x = year, y = death_rate, group = age, color = as.factor(age))) +
  geom_line() +
  labs(title = "死亡率随年份和年龄的变化",
       x = "年份",
       y = "死亡率",
       color = "年龄") +
  theme_minimal()


df = colSums(as.data.frame(Fukushima$death$total$train[,44:66]))/colSums(as.data.frame(Fukushima$expos$total$train[,44:66]))
plot(log(df), type = "l")
df


df = colSums(as.data.frame(region_Tohoku$death$total$train[,44:66]))/colSums(as.data.frame(region_Tohoku$expos$total$train[,44:66]))
df
plot(df, type = "l")

df = colSums(as.data.frame(region_Kinki$death$total$train[,44:66]))/colSums(as.data.frame(region_Kinki$expos$total$train[,44:66]))
df
plot(df, type = "l")

