library(glmnet)

set.seed(12345)

q <- 5
b0 <- 0
b <- rep(0,q)
n <- 100
N <- 1000
X <- matrix(0,n,q)
numlist <- seq(-12,-3,0.01)
numlist_lasso <- seq(-10,-3,0.1)
l <- length(numlist)
l_lasso <- length(numlist_lasso)
lambdalist <- exp(numlist)
lambdalist_lasso <- exp(numlist_lasso)
results <- vector("list", N)
results_array <- array(NA, dim = c(q+1, l, N))
ra <-  array(NA, dim = c(q+1, 1, N))
rar <- results_array
ral <- array(NA, dim = c(q+1, l_lasso, N))
rael <- results_array
rm(results_array)


b0 <- round(rnorm(1,0,1),2)
b <- round(rnorm(q,0,10),2)

#b0 <- 0
#b <- rep(10,5)


for(j in 1:N){
  
  for(i in 1:q){
    X[,i] <- rnorm(n,0,1)
  }
  Xt <- X
  
  #sds <- sqrt(colMeans(X^2) - colMeans(X)^2)
  
  #Xt <- t((t(X)-colMeans(X))/sds)
  
  eta <- b0 + Xt%*%b
  p <- plogis(eta)
  Y <- rbinom(n,1,p)
  
  gg <- glm(Y~Xt, family = "binomial")
  ggr <- glmnet(Xt,Y, family = "binomial", lambda = lambdalist, alpha = 0)
  ggl <- glmnet(Xt,Y, family = "binomial", lambda = lambdalist_lasso, alpha = 1)
  ggel <- glmnet(Xt,Y, family = "binomial", lambda = lambdalist, alpha = 0.5)
  
  ra[,,j] <- as.matrix(coef(gg))
  rar[,,j] <- as.matrix(coef(ggr))
  ral[,,j] <- as.matrix(coef(ggl))
  rael[,,j] <- as.matrix(coef(ggel))
  
  print(j/N)
}



bb <- c(b0,b)
lamlist2 <- log(ggr$lambda)
lamlist2_lasso <- log(ggl$lambda)



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

maxmse <- max(MSEr,
              MSEl,
              MSEel)

gr <- which.min(MSEr)
gl <- which.min(MSEl)
gel <- which.min(MSEel)

org <- c(v,bi2,MSE)
ridge <- c(v_r[gr], bi2_r[gr], MSEr[gr], lamlist2[gr])
lasso <- c(v_l[gl], bi2_l[gl], MSEl[gl], lamlist2_lasso[gl])
elastic <- c(v_el[gel], bi2_el[gel], MSEel[gel], lamlist2[gel])






plot(lamlist2, MSEr, type = "l", main = paste("Ridge fast. voll. Sep., n = ",n,", N = ",N, sep = ""), ylim = c(0,80),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, v_r, col = "#21918c", lwd = 2.5)
lines(lamlist2, bi2_r, col = "#fde725", lwd = 2.5)


plot(lamlist2_lasso, MSEl, type = "l", main = paste("Lasso fast. voll. Sep., n = ",n,", N = ",N, sep = ""),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2_lasso, v_l, col = "#21918c", lwd = 2.5)
lines(lamlist2_lasso, bi2_l, col = "#fde725", lwd = 2.5)


plot(lamlist2, MSEel, type = "l", main = paste("EN 0.5 fast. voll. Sep., n = ",n,", N = ",N, sep = ""),
     col = "#440154", lwd = 3, xlab = "log(Lambda)", ylab = "MSE")
lines(lamlist2, v_el, col = "#21918c", lwd = 2.5)
lines(lamlist2, bi2_el, col = "#fde725", lwd = 2.5)


plot(lamlist2, log(MSEr), type = "l", col = "#1F77B4", lwd = 2.5, ylim = c(0,10),
     main = paste("Vergleich MSE fast. voll. Sep., n = ",n,", N = ",N, sep = ""), xlab = "log(Lambda)", ylab = "log(MSE)")
lines(lamlist2_lasso, log(MSEl), col = "#FF7F0E", lwd = 2.5)
lines(lamlist2, log(MSEel), col = "#2CA02C", lwd = 2.5)




round(ridge,3)
round(lasso,3)
round(elastic,3)
round(org,3)


