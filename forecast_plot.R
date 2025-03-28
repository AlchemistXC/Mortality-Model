load("compare/region_B.RData")
load("compare/prefecture_B.RData")
load("compare/prefecture_C.RData")
mf = allresult$result$prefecture_LC_C_2.RDS$mf
# mfore_i[[i]] <- array(mfore[,i], c(35, Tfore, 8))

i = 1
j = 1

data = data.frame(m_2.5 = log(mf[[1]][i,,j]),
                  m_25 = log(mf[[2]][i,,j]),
                  Mortality_rate = log(mf[[3]][i,,j]),  
                  m_75 = log(mf[[4]][i,,j]),
                  m_97.5 = log(mf[[5]][i,,j]))


#0岁 2011  2021  2031
new_matrix <- data.frame(a0_2017 = log(mf[[3]][1, 22,]),
                         a0_2022 = log(mf[[3]][1, 32,]),
                         a0_2027 = log(mf[[3]][1, 42,]),
                         a25_2017 = log(mf[[3]][10, 22,]),
                         a25_2022 = log(mf[[3]][10, 32,]),
                         a25_2037 = log(mf[[3]][10, 42,]),
                         a55_2017 = log(mf[[3]][20, 22,]),
                         a55_2022 = log(mf[[3]][20, 32,]),
                         a55_2037 = log(mf[[3]][20, 42,]))  # Create an empty data frame with one column

new_matrix
state
real_data = array(0, c(35, 33, 46))
for (i in 1:46) {
  real_data[,,i] = get(state[i])$rate$total$raw[,44:76]
}
new_matrix <- data.frame(real_2020 = real_data[20, 31,],
                         real_2021 = real_data[20, 32,],
                         real_2022 = real_data[20, 33,],
                         fore_2020 = mf[[3]][20, 31,],
                         fore_2021 = mf[[3]][20, 32,],
                         fore_2022 = mf[[3]][20, 33,])  # Create an empty data frame with one column

new_matrix

library(sf)
library(ggplot2)
states01 <- sf::read_sf("shp//jpn.shp")
states01
library(ggplot2)

# 删除第47行
states01 <- states01[-47, ]
states01 <- cbind(states01, new_matrix)
write.csv(new_matrix, "prefecture_forecast.csv")
write.csv(new_matrix, "prefecture_covid19.csv")
# 绘制地图
ggplot(data=states01)+geom_sf()+theme_void()

states01 <- st_as_sf(states01, encoding = "UTF-8")
st_write(states01, "shp/JP_prefecture.shp")
# # 使用ggplot2绘图
# library(ggplot2)
# 
# ggplot(data = data[1:33, ]) + 
#   geom_ribbon(alpha = 0.3, aes(x = 2000:2032, ymin = m_2.5, ymax = m_97.5), fill = "#FF99FF") +
#   geom_ribbon(alpha = 0.3, aes(x = 2000:2032, ymin = m_25, ymax = m_75), fill = "#FF9966") +
#   geom_line(aes(x = 2000:2032, y = Mortality_rate))

a0 = log(mf[[3]][1, ,])
a25 = log(mf[[3]][10, ,])
a55 = log(mf[[3]][20, ,])
write.csv(a0, "a0.csv")
write.csv(a25, "a25.csv")
write.csv(a55, "a55.csv")
