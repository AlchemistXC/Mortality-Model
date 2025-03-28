state = c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki",  "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata",
          "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane",
          "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima") #, "Okinawa"

region = c("region_Hokkaido","region_Tohoku", "region_Kanto", "region_Chubu", "region_Kinki", "region_Chugoku", "region_Shikoku", "region_Kyushu")

region_list = list(c(1),2:7, 8:14, 15:23, 24:30, 31:35, 36:39, 40:46)

region_data <- function(i, nf = 10){
  year = 1947:2022
  obj <- list(type = "mortality", label = region[i])
  obj$rate <- obj$expos <- obj$death <- list()
  for (gender in 1:3) {
    obj$rate[[gender]] <- obj$expos[[gender]] <- obj$death[[gender]] <- list()
    obj$death[[gender]][[1]] <- matrix(0,nrow = 35, ncol = 76)
    obj$expos[[gender]][[1]] <- matrix(0,nrow = 35, ncol = 76)
  }
  for (j in region_list[[i]]){
    for (gender in 1:3) {
      obj$death[[gender]][[1]] <- obj$death[[gender]][[1]] + get(state[j])$death[[gender]][[1]]
      obj$expos[[gender]][[1]] <- obj$expos[[gender]][[1]] + get(state[j])$expos[[gender]][[1]]
    }
  }
  for (gender in 1:3) {
    obj$rate[[gender]][[1]] <- obj$death[[gender]][[1]] / obj$expos[[gender]][[1]]
  }
  for (gender in 1:3) {
    # train
    obj$death[[gender]][[2]] <- obj$death[[gender]][[1]][,1:(length(year)-nf)]
    obj$expos[[gender]][[2]] <- obj$expos[[gender]][[1]][,1:(length(year)-nf)]
    obj$rate[[gender]][[2]] <- obj$rate[[gender]][[1]][,1:(length(year)-nf)]
    # test
    obj$death[[gender]][[3]] <- obj$death[[gender]][[1]][,(length(year)-nf+1):length(year)]
    obj$expos[[gender]][[3]] <- obj$expos[[gender]][[1]][,(length(year)-nf+1):length(year)]
    obj$rate[[gender]][[3]] <- obj$rate[[gender]][[1]][,(length(year)-nf+1):length(year)]
    
    names(obj$death[[gender]]) = (names(obj$expos[[gender]]) = (names(obj$rate[[gender]]) <- c("raw", "train", "test")))
  }
  names(obj$expos) = (names(obj$rate) = (names(obj$death) <- tolower(c("Female", "Male", "Total" ))))
  #save
  return(obj)
}
for (i in 1:length(region)) {
  assign(region[i], region_data(i))
}
