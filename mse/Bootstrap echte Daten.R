#ALL SINGLE USE!!!!!!!!!!!!!!
library(glmnet)
heart <- read.csv("heart.xls", stringsAsFactors=TRUE)

see <- 12345

set.seed(see)

Y <- as.numeric(heart[,12])
Xdf <- heart[,-12]


mmX <- model.matrix(Y~., Xdf)[,-1]

q <- ncol(mmX)
n <- nrow(mmX)
N <- 1000
numlist <- seq(-12,3,0.1)
l <- length(numlist)
lambdalist <- exp(numlist)
results_array <- array(NA, dim = c(q+1, l, N))
ra <-  array(NA, dim = c(q+1, 1, N))
rar <- results_array
ral <- results_array
rael <- results_array
rm(results_array)

Xt <- mmX

sds <- sqrt(colMeans(Xt^2) - colMeans(Xt)^2)

Xt <- t((t(mmX)-colMeans(mmX))/sds)


gg1 <- glmnet(Xt,Y, family = "binomial", alpha = 0, lambda = 0)
org_coeffs <- coef(gg1)
bb <- org_coeffs@x




for(i in 1:N){
  

  
  samp <- sample(1:n,n,replace = T)
  Xs <- Xt[samp,]
  Ys <- Y[samp]
  
  gg <- glmnet(Xs,Ys, family = "binomial", alpha = 0, lambda = 0)
  ggr <- glmnet(Xs,Ys, family = "binomial", alpha = 0, lambda = lambdalist)
  ggl <- glmnet(Xs,Ys, family = "binomial", alpha = 1, lambda = lambdalist)
  ggel <- glmnet(Xs,Ys, family = "binomial", alpha = 0.5, lambda = lambdalist)
  
  ra[,,i] <- as.matrix(coef(gg))
  rar[,,i] <- as.matrix(coef(ggr))
  ral[,,i] <- as.matrix(coef(ggl))
  rael[,,i] <- as.matrix(coef(ggel))
  
  print(i/N * 100)
}







lamlist2 <- log(ggl$lambda)



m <- apply(ra,c(1,2),mean)
m_r <- apply(rar,c(1,2),mean)
m_l <- apply(ral,c(1,2),mean)
m_el <- apply(rael,c(1,2),mean)

bi <- m - bb
bi_r <- m_r - bb
bi_l <- m_l - bb
bi_el <- m_el - bb

bi2 <- colSums(bi^2)
bi2_r <- colSums(bi_r^2)
bi2_l <- colSums(bi_l^2)
bi2_el <- colSums(bi_el^2)

v <- colSums(apply(ra,c(1,2),var))
v_r <- colSums(apply(rar,c(1,2),var))
v_l <- colSums(apply(ral,c(1,2),var))
v_el <- colSums(apply(rael,c(1,2),var))

MSE <- v + bi2
MSEr <- v_r + bi2_r
MSEl <- v_l + bi2_l
MSEel <- v_el + bi2_el

maxmse <- max(MSEr,MSEl,MSEel)

gr <- which.min(MSEr)
gl <- which.min(MSEl)
gel <- which.min(MSEel)

org <- c(v,bi2,MSE)
ridge <- c(v_r[gr], bi2_r[gr], MSEr[gr], lamlist2[gr])
lasso <- c(v_l[gl], bi2_l[gl], MSEl[gl], lamlist2[gl])
elastic <- c(v_el[gel], bi2_el[gel], MSEel[gel], lamlist2[gel])




plot(lamlist2, MSEr, type = "l", ylim = c(0,maxmse), main = paste("Ridge heart data, n = ",n,", N = ",N, sep = ""),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, v_r, col = "#21918c", lwd = 2.5)
lines(lamlist2, bi2_r, col = "#fde725", lwd = 2.5)


plot(lamlist2, MSEl, type = "l", ylim = c(0,maxmse), main = paste("Lasso heart data, n = ",n,", N = ",N, sep = ""),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, v_l, col = "#21918c", lwd = 2.5)
lines(lamlist2, bi2_l, col = "#fde725", lwd = 2.5)


plot(lamlist2, MSEel, type = "l", ylim = c(0,maxmse), main = paste("EN 0.5 heart data, n = ",n,", N = ",N, sep = ""),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, v_el, col = "#21918c", lwd = 2.5)
lines(lamlist2, bi2_el, col = "#fde725", lwd = 2.5)


plot(lamlist2, MSEr, type = "l", ylim = c(0,maxmse), col = "#1F77B4", lwd = 2.5,
     main = paste("Vergleich MSE heart data, n = ",n,", N = ",N, sep = ""), xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, MSEl, col = "#FF7F0E", lwd = 2.5)
lines(lamlist2, MSEel, col = "#2CA02C", lwd = 2.5)




round(ridge,3)
round(lasso,3)
round(elastic,3)
round(org,3)




stop("ab hier nur noch tests")


plot(ggl)

sd <- apply(ra,c(1,2),sd)
sdr <- apply(rar,c(1,2),sd)
sdl <- apply(ral,c(1,2),sd)
sdel <- apply(rael,c(1,2),sd)

sd
sdr[,gr]
sdl[,gl]
sdel[,gel]

sdss <- data.frame(sd,sd,sdr[,gr],sdl[,gl],sdel[,gel])


bbcorrected <- bb - bi
as.numeric(bbcorrected)
bb
m_r[,gr]

plot(1:16,bb, col = "black", pch = 16)
points(1:16,as.numeric(bbcorrected), col = "red")
points(1:16,m_r[,gr], col = "yellow3", pch = 4)
points(1:16,m_l[,gl], col = "blue3", pch = 4)
points(1:16,m_el[,gel], col = "gray50", pch = 4)

dfdf <- data.frame(org = bb, corrected = as.numeric(bbcorrected), ridge = m_r[,gr], lasso = m_l[,gl], elasticnet = m_el[,gel])

upper <- dfdf + sdss*2
lower <- dfdf - sdss*2

plot(1:16,dfdf$ridge)
points(1:16,upper$ridge, col = "red")
points(1:16,lower$ridge, col = "red")

pairs(dfdf)
