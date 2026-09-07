library(glmnet)
heart <- read.csv("heart.xls", stringsAsFactors=TRUE)
Y <- as.numeric(heart[,12])
Xdf <- heart[,-12]
mmX <- model.matrix(Y~., Xdf)[,-1]

numlist <- seq(-12,3,0.1)
l <- length(numlist)
lambdalist <- exp(numlist)

n <- nrow(mmX)





ncol(mmX)
set.seed(12345)
#10f Cross-Validation
cvcvR <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0, nfolds = 10, type.measure = "deviance", lambda = lambdalist)
plot(cvcvR, sign.lambda = 1)
cvcvRm <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0, nfolds = 10, type.measure = "class", lambda = lambdalist)
plot(cvcvRm, sign.lambda = 1)

cvcvL <- cv.glmnet(mmX,Y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "deviance", lambda = lambdalist)
plot(cvcvL, sign.lambda = 1)
cvcvLm <- cv.glmnet(mmX,Y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "class", lambda = lambdalist)
plot(cvcvLm, sign.lambda = 1)

cvcvEL <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0.5, nfolds = 10, type.measure = "deviance", lambda = lambdalist)
plot(cvcvEL, sign.lambda = 1)
cvcvELm <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0.5, nfolds = 10, type.measure = "class", lambda = lambdalist)
plot(cvcvELm, sign.lambda = 1)




cvcvR
cvcvRm
cvcvL
cvcvLm
cvcvEL
cvcvELm

cvcvL
cvcvLm

round(min(cvcvR$cvm),3)
round(min(cvcvL$cvm),3)
round(min(cvcvEL$cvm),3)

round(min(cvcvRm$cvm),3)
round(min(cvcvLm$cvm),3)
round(min(cvcvELm$cvm),3)


log(cvcvR$lambda.min)
log(cvcvL$lambda.min)
log(cvcvEL$lambda.min)

log(cvcvR$lambda.1se)
log(cvcvL$lambda.1se)
log(cvcvEL$lambda.1se)

log(cvcvRm$lambda.min)
log(cvcvLm$lambda.min)
log(cvcvELm$lambda.min)

log(cvcvRm$lambda.1se)
log(cvcvLm$lambda.1se)
log(cvcvELm$lambda.1se)

#######################################################

set.seed(12345)
cvcv <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0.5, nfolds = 10, type.measure = "deviance", lambda = c(0,0.00001,0.00002,1))
cvcv$cvm

set.seed(12345)
cvcvc <- cv.glmnet(mmX,Y, family = "binomial", alpha = 0.5, nfolds = 10, type.measure = "class", lambda = c(0,0.00001,0.00002,1))
cvcvc$cvm

