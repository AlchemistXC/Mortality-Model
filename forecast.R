load("compare/region_B.RData")

mu = allresult$result$region_LC_C.RDS$mfore
state = c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki",  "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata",
          "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane",
          "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
yearf <- 2023:2032
length(yearf)

data = data.frame(Time = 1:10,
                  m_2.5 = log(mu[[1]][3,,3]),
                  m_25 = log(mu[[2]][3,,3]),
                  Moratility_rate = log(mu[[3]][3,,3]),
                  m_75 = log(mu[[4]][3,,3]),
                  m_97.5 = log(mu[[5]][3,,3]))
ggplot() + 
  geom_ribbon(data = data, alpha = 0.3, aes(x = Time, ymin = m_2.5, ymax = m_97.5), fill = "#FF99FF") +
  geom_ribbon(data = data, alpha = 0.3, aes(x = Time, ymin = m_25, ymax = m_75), fill = "#FF9966") +
  geom_line(data = data, aes(x = Time, y = Moratility_rate)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
