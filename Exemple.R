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


#local.dir = "~/Documents/BLN"
local.dir = getwd()

sourceCpp(file.path(local.dir,"src/optim_rank_cov_BLN.cpp"))

# Simulation d'un jeu de données simple

n <- 100
p <- 10
d <- 3
q <- 5

X <- cbind(c(rep(1, n*p)),matrix(rnorm(n*p*d), nrow = n*p))
B <- c(-2, 1, 10, 2)

W <- matrix(rnorm(n*q), nrow = n)
C <- matrix(rnorm(p*q), nrow = p)/sqrt(q)


XB <- VectorToMatrix(X%*%B, n, p)
Z <- XB + W%*%t(C)    


Prob <- plogis(Z)

Y <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)

plot(Z, Y)
lines(sort(Z), plogis(sort(Z)))
boxplot(Z ~ Y)


# Initialisation avec une régression logistique

vecY <- MatrixToVector(Y)

fit <- glm(vec(Y) ~ -1 + X, family = "binomial")

B.init <- as.matrix(fit$coefficients)
#B.init <- matrix(c(0, 0, 5, 1), ncol =1)
res.vec <- fit$residuals
res.full <- ifelse(is.na(vecY), 0, res.vec)
res.mat <- VectorToMatrix(res.full, n, p)

svdM <- svd(res.mat, nu = q, nv = p)

C.init <- svdM$v[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q)/sqrt(n)
M.init  <- svdM$u[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q) %*% t(svdM$v[1:q, 1:q, drop = FALSE])
S.init <-  matrix(0.1, n, q)

Mu.init <- VectorToMatrix(X%*%B.init, n, p)

# C.init <- matrix(1, p, q)
# M.init <- matrix(0, n, q)


# Application du jeu de données


params <- list(B = B.init, C = C.init, M = M.init, S = S.init)

data <- list(Y = Y,
             R = matrix(1, n, p),
             X = X)

config <- PLNPCA_param()$config_optim

# sourceCpp(file.path(local.dir,"src/optim_rank_cov_BLN.cpp"))

out <- nlopt_optimize_rank_BLN(data, params, config)


plot(B, out$B,xlim = c(-2,10),ylim = c(-6,12)); abline(0,1) ; points(B, B.init, col = "red")
plot(C%*%t(C), out$C%*%t(out$C)) ; abline(0,1)

plot(XB, Mu.init); abline(0,1)
plot(Z, VectorToMatrix(X%*%B.init, n, p) + M.init%*%t(C.init)); abline(0,1)
plot(Z, VectorToMatrix(X%*%B.init, n, p)); abline(0,1)
plot(Z, VectorToMatrix(X%*%out$B, n, p) + out$M%*%t(out$C)); abline(0,1)
Z.pred <- VectorToMatrix(X%*%out$B, n, p) + out$M%*%t(out$C)
mean(Z.pred)




plot(C%*%t(C), out$C%*%t(out$C)); abline(0,1)


#-------------------------------------------------------------------------------
# Z = Mu.init + M.init %*% t(C.init)
# 
# expA = exp(A)
# expmA = exp(-A)
# 
# D <- 1 / (1 + expmA)
# E <- (1 - expmA)/ (4*hadamard.prod(A, (1 + expmA)))
# f <- Y + D + 2*hadamard.prod(E, (Mu.init + M.init %*% t(C.init) - A))
# G1 <- log(1 + expA) - hadamard.prod(E, A)
# G2 <- log(1 + expmA) + A - hadamard.prod(E, A)
#   
# el1 <- hadamard.prod((Y + D), A1)
# el2 <- hadamard.prod(E, A2^2 + (A1 - A)^2)
# obj <- sum()
  








