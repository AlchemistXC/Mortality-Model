state = c("Japan","Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki",  "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata",
          "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane",
          "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima") # "Okinawa"
# region = c("ind_R1", "ind_R2", "ind_R3", "ind_R4", "ind_R5", "ind_R6", "ind_R7", "ind_R8")

read.jpn1_3 <- function (region = "01",  label, fore = 10)
{
  filename_deaths = paste("data/",region,"_deaths",".txt",sep = "")
  if(!file.exists(filename_deaths)){
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Deaths_1x1.txt", sep = "")
    utils::download.file(url=path,destfile=filename_deaths)
  }
  dt <- try(read.table(filename_deaths, skip = 2, header = TRUE, na.strings = "\n"), TRUE)
  if(class(dt) == "try-error")
    stop("Connection error at www.mortality.org. Please check username, password and country label.")
  filename_exposures = paste("data/",region,"_exposures",".txt",sep = "")
  if(!file.exists(filename_exposures)){
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Exposures_1x1.txt", sep = "")
    utils::download.file(url=path,destfile=filename_exposures)
  }
  pop <- try(read.table(filename_exposures, skip = 2, header = TRUE, na.strings = "."),
             TRUE)
  if (class(pop) == "try-error")
    stop("Exposures file not found at www.mortality.org")
  obj <- list(type = "mortality", label = label)
  n <- length(unique(dt[, 1]))
  m <- length(unique(dt[, 2]))
  age <- c("0")
  for (i in 2:35) {
    if(i != 35){age[i] <- paste(3*(i-2)+1,"-",3*(i-1),sep = "")}
     else{age[i] = "100+"}
  }
  year <- 1947:2022
  mnames <- names(dt)[-c(1, 2)]
  m2names <- c("raw", "train", "test")
  n.mort <- length(mnames)
  obj$rate <- obj$expos <- obj$death <- list()
  for (i in 1:n.mort) {
    obj$rate[[i]] <- obj$expos[[i]] <- obj$death[[i]] <- list()
    # raw
    a <- matrix(dt[, i + 2], nrow = m, ncol = n)
    b <- matrix(pop[, i + 2], nrow = m, ncol = n)
    foodeath <- matrix(0,nrow = length(age), ncol = n)
    fooexpos <- matrix(0,nrow = length(age), ncol = n)
    obj$death[[i]][[1]] <- matrix(0,nrow = length(age), ncol = length(year))
    obj$expos[[i]][[1]] <- matrix(0,nrow = length(age), ncol = length(year))
    for (j in 1:length(age)) {
      if(j == 1){
        obj$death[[i]][[1]][j,] <- a[1,]
        obj$expos[[i]][[1]][j,] <- b[1,]
      }
      else{
        if(j < length(age)){
          obj$death[[i]][[1]][j,] <- colSums(a[(3*(j-2)+1):(3*(j-1)),])
          obj$expos[[i]][[1]][j,] <- colSums(b[(3*(j-2)+1):(3*(j-1)),])
        }
        else{
          obj$death[[i]][[1]][j,] <- colSums(a[(3*(j-2)+1):m,])
          obj$expos[[i]][[1]][j,] <- colSums(b[(3*(j-2)+1):m,])
        }
      }

    }
    obj$death[[i]][[1]] <- floor(obj$death[[i]][[1]])
    obj$expos[[i]][[1]] <- apply(obj$expos[[i]][[1]], c(1, 2), function(x) {
      if (x == 0) {
        return(1)
      } else {
        return(ceiling(x))
      }
    })
    #obj$expos[[i]][[1]] <- floor(obj$expos[[i]][[1]])
    obj$rate[[i]][[1]] <- obj$death[[i]][[1]] / obj$expos[[i]][[1]]
    dimnames(obj$death[[i]][[1]]) <- dimnames(obj$expos[[i]][[1]]) <- dimnames(obj$rate[[i]][[1]]) <- list(age,year)

    # train
    nf <- floor(fore)
    obj$death[[i]][[2]] <- obj$death[[i]][[1]][,1:(length(year)-nf)]
    obj$expos[[i]][[2]] <- obj$expos[[i]][[1]][,1:(length(year)-nf)]
    obj$rate[[i]][[2]] <- obj$rate[[i]][[1]][,1:(length(year)-nf)]
    # test
    obj$death[[i]][[3]] <- obj$death[[i]][[1]][,(length(year)-nf+1):length(year)]
    obj$expos[[i]][[3]] <- obj$expos[[i]][[1]][,(length(year)-nf+1):length(year)]
    obj$rate[[i]][[3]] <- obj$rate[[i]][[1]][,(length(year)-nf+1):length(year)]
    #name("raw", "train", "test")
    names(obj$death[[i]]) = (names(obj$expos[[i]]) = (names(obj$rate[[i]]) <- m2names))
  }
  names(obj$expos) = (names(obj$rate) = (names(obj$death) <- tolower(mnames)))
  return(obj)
}

# read.jpn1_3 <- function (region = "01",  label, fore = 10)
# {
#   filename_deaths = paste("data/",region,"_deaths",".txt",sep = "")
#   if(!file.exists(filename_deaths)){
#     path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Deaths_1x1.txt", sep = "")
#     utils::download.file(url=path,destfile=filename_deaths)
#   }
#   dt <- try(read.table(filename_deaths, skip = 2, header = TRUE, na.strings = "\n"), TRUE)
#   if(class(dt) == "try-error")
#     stop("Connection error at www.mortality.org. Please check username, password and country label.")
#   filename_exposures = paste("data/",region,"_exposures",".txt",sep = "")
#   if(!file.exists(filename_exposures)){
#     path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Exposures_1x1.txt", sep = "")
#     utils::download.file(url=path,destfile=filename_exposures)
#   }
#   pop <- try(read.table(filename_exposures, skip = 2, header = TRUE, na.strings = "."),
#              TRUE)
#   if (class(pop) == "try-error")
#     stop("Exposures file not found at www.mortality.org")
#   obj <- list(type = "mortality", label = label)
#   n <- length(unique(dt[, 1]))
#   m <- length(unique(dt[, 2]))
#   age <- c("0")
#   for (i in 2:35) {
#     age[i] <- paste(3*(i-2)+1,"-",3*(i-1),sep = "")
#   }
#   year <- 1947:2022
#   mnames <- names(dt)[-c(1, 2)]
#   m2names <- c("raw", "train", "test")
#   n.mort <- length(mnames)
#   obj$rate <- obj$expos <- obj$death <- list()
#   for (i in 1:n.mort) {
#     obj$rate[[i]] <- obj$expos[[i]] <- obj$death[[i]] <- list()
#     # raw
#     a <- matrix(dt[, i + 2], nrow = m, ncol = n)
#     b <- matrix(pop[, i + 2], nrow = m, ncol = n)
#     foodeath <- matrix(0,nrow = length(age), ncol = n)
#     fooexpos <- matrix(0,nrow = length(age), ncol = n)
#     obj$death[[i]][[1]] <- matrix(0,nrow = length(age), ncol = length(year))
#     obj$expos[[i]][[1]] <- matrix(0,nrow = length(age), ncol = length(year))
#     for (j in 1:length(age)) {
#       if(j == 1){
#         obj$death[[i]][[1]][j,] <- a[1,]
#         obj$expos[[i]][[1]][j,] <- b[1,]
#       }
#       else{
#         obj$death[[i]][[1]][j,] <- colSums(a[(3*(j-2)+1):(3*(j-1)),])
#         obj$expos[[i]][[1]][j,] <- colSums(b[(3*(j-2)+1):(3*(j-1)),])
#       }
#       
#     }
#     obj$death[[i]][[1]] <- floor(obj$death[[i]][[1]])
#     obj$expos[[i]][[1]] <- apply(obj$expos[[i]][[1]], c(1, 2), function(x) {
#       if (x == 0) {
#         return(1)
#       } else {
#         return(ceiling(x))
#       }
#     })
#     #obj$expos[[i]][[1]] <- floor(obj$expos[[i]][[1]])
#     obj$rate[[i]][[1]] <- obj$death[[i]][[1]] / obj$expos[[i]][[1]]
#     dimnames(obj$death[[i]][[1]]) <- dimnames(obj$expos[[i]][[1]]) <- dimnames(obj$rate[[i]][[1]]) <- list(age,year)
#     
#     # train
#     nf <- floor(fore)
#     obj$death[[i]][[2]] <- obj$death[[i]][[1]][,1:(length(year)-nf)]
#     obj$expos[[i]][[2]] <- obj$expos[[i]][[1]][,1:(length(year)-nf)]
#     obj$rate[[i]][[2]] <- obj$rate[[i]][[1]][,1:(length(year)-nf)]
#     # test
#     obj$death[[i]][[3]] <- obj$death[[i]][[1]][,(length(year)-nf+1):length(year)]
#     obj$expos[[i]][[3]] <- obj$expos[[i]][[1]][,(length(year)-nf+1):length(year)]
#     obj$rate[[i]][[3]] <- obj$rate[[i]][[1]][,(length(year)-nf+1):length(year)]
#     #name("raw", "train", "test")
#     names(obj$death[[i]]) = (names(obj$expos[[i]]) = (names(obj$rate[[i]]) <- m2names))
#   }
#   names(obj$expos) = (names(obj$rate) = (names(obj$death) <- tolower(mnames)))
#   return(obj)
# }

for (i in 1:length(state)) {
  if(i < 11){index = paste("0", as.character(i-1), sep = "")}
  else{index = as.character(i-1)}
  assign(state[i], read.jpn1_3(index, state[i]))
  print(i)
}

