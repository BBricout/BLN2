library(Rcpp)
library(RcppArmadillo)
library(nloptr)

###################################""
# From matrix to vector
##################################
MatrixToVector <- function(matrix) {
  n <- nrow(matrix)
  p <- ncol(matrix)
  # Check if the input matrix is not empty
  if (n == 0 || p == 0) {
    stop("Input matrix is empty")
  }
  # Reshape the matrix into a column vector
  vectorized <- as.vector(matrix)
  return(vectorized)
}

####################################
# From vector to matrix
######################################
VectorToMatrix <- function(vector, n, p) {
  # Create a matrix with the specified dimensions
  matrix <- matrix(vector, nrow = n, ncol = p)
  return(matrix)
}


###############################"
# From vector to parametersList
##############################

VectorToParamsList <- function(params, list_dim) {
  
  
  d <- list_dim$d
  n <- list_dim$n
  p <- list_dim$p
  q <- list_dim$q
  
  B <- matrix(params[1:d],ncol=1)
  
  C <- VectorToMatrix(params[(d + 1):(d + p*q)],p,q)
  M <- VectorToMatrix(params[(d+ p*q  + 1 ):(d + p*q + n*q)],n,q)
  S <- VectorToMatrix(params[(d+ p*q  + n*q +1 ):(d + p*q + 2*n*q)],n,q)
  return(list(B=B,C=C,M=M,S=S))
}

###############################"
# From vector to parametersList
##############################

ParamsListToVector <- function(paramsList, list_dim) {
  
  d <- list_dim$d
  n <- list_dim$n
  p <- list_dim$p
  q <- list_dim$q
  
  params <-  rep(0,d + p*q + 2*n*q)
  params[1:d] <- c(paramsList$B)
  params[(d + 1):(d + p*q)] <- MatrixToVector(paramsList$C)
  params[(d+ p*q  + 1 ):(d + p*q + n*q)] <- MatrixToVector(paramsList$M)
  params[(d+ p*q  + n*q +1 ):(d + p*q + 2*n*q)] <- MatrixToVector(paramsList$S)
  return(params)
}


#################################################
objective<- function(params, data) {
  
  Y <- data$Y
  X <- data$X
  R <- data$R
  
  list_dim <- list(n = 500, p = 30, d = 6, q = 3)
  
  paramslist <- VectorToParamsList(params, list_dim)
  
  B <- paramslist$B
  C <- paramslist$C
  M <- paramslist$M
  S <- paramslist$S
  
  

  #---------------- Evaluation objectif
  XB <- data$X %*% B
  Mu <- VectorToMatrix(XB, n,p)
  
  A <- sqrt((Mu + M %*% t(C))^2 + (S %*% t(C*C))^2)
  expmA <- exp(-A)
  D <- 1 / (1 + expmA)
  E <- (1-expmA) / (4 * A * (1+expmA))
  G <- log(1 + expmA) + A  - E * A
  
  objective <- sum(R * ((Y - D) * (Mu + M %*% t(C))) - E * ((S %*% t(C*C)) + (Mu + M %*% t(C) - A)^2) - G) - 0.5 * sum((M^2 + S - log(S)))
  objective <- -1*objective 
  return(objective)
}


nlopt_optimize_rank_BLN <- function(data, params, eval_f, list_dim) {
  init_params <- ParamsListToVector(params, list_dim)
  result <- nloptr(x0 = init_params, eval_f = eval_f, 
                   opts = list("algorithm" = "NLOPT_GN_ISRES", "xtol_rel" = 1.0e-8),
                   data = data)
  return(result)
}

list_dim <- list(n = 500, p = 30, d = 6, q = 3)

res <- nlopt_optimize_rank_BLN(data, params, objective, list_dim)

VectorToParamsList(res$solution, list_dim)$C

































################################################



# # Rank-constrained covariance
# nlopt_optimize_rank_BLN <- function(data, params, eval_f) {
#   # Conversion from R, prepare optimization
#   Y <- as.matrix(data$Y)  # responses (n,p)
#   R <- as.matrix(data$R)  # missing data (n,p)
#   X <- as.matrix(data$X)  # covariates (np,d)
#   
#   init_S <- params$S
#   init_B <- params$B
#   init_M <- params$M
#   init_C <- params$C
#   
#   data <- list(Y = Y, X = X, R = R, B = init_B, C = init_C, M = init_M, S = init_S)
#   
#   
#   # Optimize
#   
#   result <- nloptr(x0 = data, eval_f = eval_f, 
#                    opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8))
#   
#   return(result)
#   
# }
objectiveB<- function(B,list_dim,data,params.fixed) {
  
  n <- list_dim$n
  p <- list_dim$p
  vecY <- MatrixToVector(data$Y)
  vecR <- MatrixToVector(data$R)
  
  #-------------- transfo paramsVector to List
  M <-  params.fixed$M
  C <-  params.fixed$C
  S <-  params.fixed$S
  #---------------- Evaluation objectif
  S2 <- S * S
  
  XB <- data$X %*% B
  Mu <- VectorToMatrix(XB, n,p)
  
  A <- 1+exp(Mu)
  
  
  #A <- sqrt((Mu + M %*% t(C))^2 + (S2 %*% t(C*C))^2)
  #expmA <- exp(-A)
  #D <- 1 / (1 + expmA)
  #E <- (1-expmA) / (4 * A * (1+expmA))
  #G <- log(1 + expmA) + A  - E * A
  
  #objective <- sum(data$R * dbinom(data$Y,prob =1/A,size = 1, log=  TRUE))
  objective <- sum(data$R * (data$Y*Mu - log(A)))
  #objective <- sum(data$R * ((data$Y + D) * (Mu + M %*% t(C))) + E * ((S2 %*% t(C*C))^2 + (Mu + M %*% t(C) - A)^2) + G) #- 0.5 * sum((M^2 + S2 + log(S2)))
  objective <- -1*objective 
  return(objective)
}


