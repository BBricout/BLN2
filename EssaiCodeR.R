rm(list=ls())
#----------Librairies-----------------

library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missForest)
library(purrr)
library(PLNmodels)


# local.dir = "~/Documents/BLN"
local.dir = getwd()

#sourceCpp(file.path(local.dir,"src/optim_rank_cov_BLN.cpp"))
source(file="codeRbyChatGPT.R")
source(file='BLNAddFunctions.R')
#---------------------------------------------------
# Simulation d'un jeu de données simple
#----------------------------------------------------
n <- 100
p <- 10
d <- 2
q <- 5

X <- matrix(rnorm(n*p*d), nrow = n*p)
X[,1] <- 1
B <- c(-2, 1) #, 10, 2)

W <- matrix(rnorm(n*q), nrow = n)
C <- 0*matrix(rnorm(p*q), nrow = p)/sqrt(q)

params.true = list(B=B,C=C)

XB <- VectorToMatrix(X%*%B, n, p)
Z <- XB + W%*%t(C)    


Prob <- plogis(Z)

Y <- matrix(rbinom(n*p, prob = Prob, size = 1), nrow = n)
Y.na <- prodNA(Y, 0.1)

plot(Z, Y)
lines(sort(Z), plogis(sort(Z)))
boxplot(Z ~ Y)


####################################################"
# Initialisation avec une régression logistique
###################################################
vecY <- MatrixToVector(Y)
fit <- glm(vecY ~ -1 + X, family = "binomial")

#B.init <- as.matrix(fit$coefficients)
B.init <- matrix(1,nrow= d, ncol =1)

res.vec <- fit$residuals
res.full <- ifelse(is.na(vecY), 0, res.vec)
res.mat <- VectorToMatrix(res.full, n, p)
svdM <- svd(res.mat, nu = q, nv = p)

# C.init <-  svdM$v[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q)/sqrt(n)
# M.init  <- svdM$u[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q) %*% t(svdM$v[1:q, 1:q, drop = FALSE])
# S.init <-  matrix(runif(n*q,0.02,0.2), n, q)

C.init <- matrix(0, p, q)
M.init <- matrix(0, n, q)
S.init <- matrix(0, n, q)

# C.init <- matrix(1, p, q)
# M.init <- matrix(0, n, q)


list_dim <- list(n=n,d=d,p=p,q=q)

params.init <- ParamsListToVector( list(B = B.init, C=C,S = S.init,M  = M.init),list_dim)

params.fixed = list(C=C,M = M.init,S=S.init)

#######################################"
# Inférence par VEM
#######################################
data <- list(Y = Y,
             R = ifelse(is.na(Y), 0, 1),
             X = X)

res.optimB <- optim(B.init,objectiveB,data = data,list_dim = list_dim,params.fixed = params.fixed)
res.optimB$par
params.true$B
fit$coefficients


