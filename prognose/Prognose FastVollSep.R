library(glmnet)

set.seed(12345)

q <- 5
b0 <- 0
b <- rep(0,q)
n <- 100
N <- 100
X <- matrix(0,n,q)
numlist <- seq(-12,3,0.1)
#l <- length(numlist)
lambdalist <- exp(numlist)
results <- array(NA, dim = c(4, N))
rar <- results
ral <- rar
rael <- rar

min_measure <- array(NA,dim=c(8,N))


b0 <- round(rnorm(1,0,1),2)
b <- round(rnorm(q,0,10),2)


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
  
  cvR <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 0)
  cvRm <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 0, type.measure = "class")
  rar[1,j] <- cvR$lambda.min
  rar[2,j] <- cvR$lambda.1se
  rar[3,j] <- cvRm$lambda.min
  rar[4,j] <- cvRm$lambda.1se
  
  cvL <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 1)
  cvLm <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 1, type.measure = "class")
  ral[1,j] <- cvL$lambda.min
  ral[2,j] <- cvL$lambda.1se
  ral[3,j] <- cvLm$lambda.min
  ral[4,j] <- cvLm$lambda.1se
  
  cvEL <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 0.5)
  cvELm <- cv.glmnet(Xt,Y, lambda = lambdalist, family = "binomial", alpha = 0.5, type.measure = "class")
  rael[1,j] <- cvEL$lambda.min
  rael[2,j] <- cvEL$lambda.1se
  rael[3,j] <- cvELm$lambda.min
  rael[4,j] <- cvELm$lambda.1se
  
  
  
  
  cvcheck <- cv.glmnet(Xt,Y, lambda = c(0,0.0000001), family = "binomial", alpha = 0.5)
  cvcheckm <- cv.glmnet(Xt,Y, lambda = c(0,0.0000001), family = "binomial", alpha = 0.5, type.measure = "class")
  
  
  min_measure[1,j] <- min(cvR$cvm)
  min_measure[2,j] <- min(cvL$cvm)
  min_measure[3,j] <- min(cvEL$cvm)
  
  min_measure[4,j] <- min(cvcheck$cvm[2])
  
  
  min_measure[5,j] <- min(cvRm$cvm)
  min_measure[6,j] <- min(cvLm$cvm)
  min_measure[7,j] <- min(cvELm$cvm)
  
  min_measure[8,j] <- min(cvcheckm$cvm[2])
  
  
  
  print(j/N)
}

med <- function(x){
  le <- length(x)
  m <- round((le+1)/2,0)
  xm <- sort(x)
  return(xm[m])
}

names <- c("dev.min", "dev.1se", "misc.min", "misc.1se")

ri <- as.data.frame(t(log(rar)))
names(ri) <- names

la <- as.data.frame(t(log(ral)))
names(la) <- names

el <- as.data.frame(t(log(rael)))
names(el) <- names

#boxplot(ri, main = "Optimale Prognose Ridge, n = 100, N = 100, fast. voll. Sep.", ylab = "log(lambda)")
#boxplot(la, main = "Optimale Prognose Lasso, n = 100, N = 100, fast. voll. Sep.", ylab = "log(lambda)")
#boxplot(el, main = "Optimale Prognose EN 0.5, n = 100, N = 100, fast. voll. Sep.", ylab = "log(lambda)")

plot(cvR, sign.lambda = 1)
plot(cvL, sign.lambda = 1)
plot(cvEL, sign.lambda = 1)

plot(cvRm, sign.lambda = 1)
plot(cvLm, sign.lambda = 1)
plot(cvELm, sign.lambda = 1)

apply(rar,1,med)
log(apply(rar,1,med))

apply(ral,1,med)
log(apply(ral,1,med))

apply(rael,1,med)
log(apply(rael,1,med))


round(matrix(apply(min_measure,1,med),4,2),3)

