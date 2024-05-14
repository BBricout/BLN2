library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missForest)

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_cov_BLN.cpp"))

#------------Simulation d'un jeu de données-------------------

n <- 500
p <- 30
d <- 5
q <- 3

X <- cbind(c(rep(1, n*p)),matrix(rnorm(n*p*d), nrow = n*p))
B <- c(-2, 1, 3, 2, 5, -3)
sum(B^2)

rho <- 1

W <- matrix(rnorm(n*q), nrow = n)
C <- matrix(rnorm(p*q,0,1), nrow = p)/sqrt(q)


XB <- VectorToMatrix(X%*%B, n, p)
Z <- XB + W%*%t(C)    


Prob <- plogis(Z)

Y <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)

data <- list(Y = Y, R = matrix(1, n, p), X = X)

var(X%*%B)
var(MatrixToVector(W%*%t(C)))

plot(Z, Prob)
boxplot(Prob~ Y)

Y.na10 <- prodNA(Y, 0.1)
Y.na30 <- prodNA(Y, 0.3)
Y.na50 <- prodNA(Y, 0.5)



O <- VectorToMatrix(X%*%B, n, p)
j=1
mat <- t(sapply(1:p, function(j)
       glm(Y[,j]~-1 + W + offset(O[,j]), family= "binomial")$coefficients))

plot(C, mat);abline(0,1)














