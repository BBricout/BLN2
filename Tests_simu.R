rm(list=ls())

library(PLNmodels)

local.dir = getwd()
source(file.path(local.dir,"Simu_data.R"))

config <- PLNPCA_param()$config_optim



# Initialisation

Init_BLN <- function(Y, X){
  fit <- glm(vec(Y) ~ -1 + X, family = "binomial", na.action = na.exclude)
  
  res <- residuals(fit, type = "response")

  ind.res <- as.numeric(names(res))
  indices <- c(1:(n*p))
  res.full <- ifelse(indices %in% ind.res, res, 0)
  
  
  B.init <- as.matrix(fit$coefficients)
  # B.init <- as.matrix(B)
  # res.mat <- VectorToMatrix(res.full, n, p)
  
  # svdM <- svd(res.mat, nu = q, nv = p)
  
  #C.init <- svdM$v[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q)/sqrt(n)
  C.init <- matrix(rnorm(p*q), nrow = p)
  # C.init <- C
  M.init <- matrix(0.1, n, q)
  # M.init <- W
  # M.init  <- svdM$u[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q) %*% t(svdM$v[1:q, 1:q, drop = FALSE])
  S.init <-  matrix(1, n, q)
  
  result <- list(B = B.init, C = C.init, M = M.init, S = S.init)
  
  return(result)
}

# Paramètres initiaux

params <- Init_BLN(Y, X)
params.na10 <- Init_BLN(Y.na10,X)
params.na30 <- Init_BLN(Y.na30,X)
params.na50 <- Init_BLN(Y.na50, X)

# Données

data <- list(Y = Y, R = matrix(1, n, p), X = X)
data.na10 <- list(Y = Y, R = ifelse(is.na(Y.na10), 0, 1), X = X)
data.na30 <- list(Y = Y, R = ifelse(is.na(Y.na30), 0, 1), X = X)
data.na50 <- list(Y = Y, R = ifelse(is.na(Y.na50), 0, 1), X = X)

# Ajustement

out <- nlopt_optimize_rank_BLN(data, params, config)
out.na10 <- nlopt_optimize_rank_BLN(data.na10, params.na10, config)
out.na30 <- nlopt_optimize_rank_BLN(data.na30, params.na30, config)
out.na50 <- nlopt_optimize_rank_BLN(data.na50, params.na50, config)

# Visualisation des regresseurs

par(mfrow = c(2, 2))  # 2 lignes, 2 colonnes

plot(B, out$B, main = "B en données complètes") ; abline(0,1)
plot(B, out.na10$B, main = "B pour 10% de NA") ; abline(0,1)
plot(B, out.na30$B, main = "B pour 30% de NA") ; abline(0,1)
plot(B, out.na50$B, main = "B pour 50% de NA") ; abline(0,1)


# Visualisation de l'ACP

par(mfrow = c(2,2))

plot(C%*%t(C), out$C%*%t(out$C), main = "C*t(C) données complètes") ; abline(0,1)
plot(C%*%t(C), out.na10$C%*%t(out.na10$C), main = "C*t(C) 10% de NA") ; abline(0,1)
plot(C%*%t(C), out.na30$C%*%t(out.na30$C), main = "C*t(C) 30% de NA") ; abline(0,1)
plot(C%*%t(C), out.na50$C%*%t(out.na50$C), main = "C*t(C) 50% de NA") ; abline(0,1)


# Visualisation de Z

par(mfrow = c(2,2))

plot(Z, VectorToMatrix(X%*%out$B, n, p) + out$M%*%t(out$C), main = "Z données complètes"); abline(0,1)
plot(Z, VectorToMatrix(X%*%out.na10$B, n, p) + out.na10$M%*%t(out.na10$C), main = "Z 10% de NA"); abline(0,1)
plot(Z, VectorToMatrix(X%*%out.na30$B, n, p) + out.na30$M%*%t(out.na30$C), main = "Z 30% de NA"); abline(0,1)
plot(Z, VectorToMatrix(X%*%out.na50$B, n, p) + out.na50$M%*%t(out.na30$C), main = "Z 50% de NA"); abline(0,1)

par(mfrow = c(2,2))

plot(VectorToMatrix(X%*%B, n, p), VectorToMatrix(X%*%out$B, n, p), main = "Complet"); abline(0,1)
plot(VectorToMatrix(X%*%B, n, p), VectorToMatrix(X%*%out.na10$B, n, p), main = "10%"); abline(0,1)
plot(VectorToMatrix(X%*%B, n, p), VectorToMatrix(X%*%out.na30$B, n, p), main = "30%"); abline(0,1)
plot(VectorToMatrix(X%*%B, n, p), VectorToMatrix(X%*%out.na50$B, n, p), main = "50%"); abline(0,1)

par(mfrow = c(1,1))

plot(W, out$M, main = "Complet") ; abline(0,1)
curve(dnorm(x, out$M[1,1], sqrt(out$S[1,1])), -10, 10, col = "red")
curve(dnorm(x, params$M[1,1], sqrt(params$S[1,1])), add = TRUE, col = "blue")
curve(dnorm(x), add = TRUE, col = "green"); abline(v=W[1,1])

plot(params$M, W)
U = as.vector(params$M[params$M < 20])
V = as.vector(W[params$M < 20]) 
plot(U, V)
var(X%*%B)
var(MatrixToVector(W%*%t(C)))
