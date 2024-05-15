library(combinat)

# Initialisation of the inference

Init_BLN <- function(data){

  n <- nrow(data$Y)
  p <- ncol(data$Y)
  
  
  Y <- MatrixToVector(data$Y)
  names(Y)<- c(1:(n*p))
  rownames(data$X) <- c(1:(n*p))
  
  w <- which(MatrixToVector(data$R)==0)
  if(length(w)>1){
    partX <- data$X[-w,]
    partY <- Y[-w]
  }else{ 
    partX <- data$X
    partY <- Y
    }
  
  fit <- glm(partY ~ -1 + partX, family = "binomial", na.action = na.exclude)
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
  M.init <- matrix(0, n, q)
  # M.init <- W
  # M.init  <- svdM$u[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q) %*% t(svdM$v[1:q, 1:q, drop = FALSE])
  S.init <-  matrix(1, n, q)
  
  result <- list(B = B.init, C = C.init, M = M.init, Sd = sqrt(S.init))
  
  return(result)
}

###################################### 
#######""
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

VectorToparamsList <- function(params.vec, list_dim,estim='all',params.list.true=NULL) {
  
  
  d <- list_dim$d
  n <- list_dim$n
  p <- list_dim$p
  q <- list_dim$q
  
  if('all' %in% estim){
    B <- params.vec[1:d]
    C <- VectorToMatrix(params.vec[(d + 1):(d + p*q)],p,q)
    M <- VectorToMatrix(params.vec[(d+ p*q  + 1 ):(d + p*q + n*q)],n,q)
    Sd <- VectorToMatrix(params.vec[(d+ p*q  + n*q +1 ):(d + p*q + 2*n*q)],n,q)
    res <- list(B=B,C=C,M=M,Sd = Sd)
  }
  if(estim == 'B'){
    res <- params.list.true
    res$B <-  params.vec[1:d]
  }
  if(estim == 'CSM'){
    C <- VectorToMatrix(params.vec[( 1):( p*q)],p,q)
    M <- VectorToMatrix(params.vec[( p*q  + 1 ):( p*q + n*q)],n,q)
    Sd <- VectorToMatrix(params.vec[( p*q  + n*q +1 ):( p*q + 2*n*q)],n,q)
    res <- list(C=C,M=M,Sd = Sd)
    res$B <- params.list.true$B
  }
  if(estim == 'C'){
    res <- params.list.true
    res$C <-  VectorToMatrix(params.vec[( 1):( p*q)],p,q)
  }
  
  return(res)
}

###############################"
# From vector to parametersList
##############################

paramListToVector <- function(params.list,list_dim,estim='all') {
  
  d <- list_dim$d
  n <- list_dim$n
  p <- list_dim$p
  q <- list_dim$q
  
  if(estim == 'all'){
    params.vec <-  rep(0,d + p*q + 2*n*q)
    params.vec[1:d] <- c(params.list$B)
    params.vec[(d + 1):(d + p*q)] <- MatrixToVector(params.list$C)
    params.vec[(d+ p*q  + 1 ):(d + p*q + n*q)] <- MatrixToVector(params.list$M)
    params.vec[(d+ p*q  + n*q +1 ):(d + p*q + 2*n*q)] <- MatrixToVector(params.list$Sd)
  }
  if(estim == 'B'){
    params.vec <- c(params.list$B)
  }
  if(estim %in% sapply(permn(c('C','S','M')),function(l){paste(l,collapse="")}))
{
    params.vec <-  rep(0,p*q + 2*n*q)
    params.vec[1:  (p*q)] <- MatrixToVector(params.list$C)
    params.vec[( p*q  + 1 ):( p*q + n*q)] <- MatrixToVector(params.list$M)
    params.vec[( p*q  + n*q +1 ):(  p*q + 2*n*q)] <- MatrixToVector(params.list$Sd)
  }
  if(estim == 'C'){
    params.vec <-  MatrixToVector(params.list$C)
  }
  return(params.vec)
}

#######################################
# Objective 
#######################################
objective<- function(params.vec, data,list_dim,estim='all',params.list.true=NULL){
  
  if(is.null(estim)){estim="all"}
  
  #
  p <- list_dim$p 
  q <- list_dim$q
  n <- list_dim$n
  d <- list_dim$d
  
  
  
  
  #print(length(params.vec))
  Y <- data$Y
  X <- data$X
  R <- data$R
  

  params.list <- VectorToparamsList(params.vec, list_dim,estim=estim,params.list.true=params.list.true)
  
  B <- params.list$B
  C <- params.list$C
  M <- params.list$M
  Sd <- params.list$Sd
  S <- params.list$Sd^2
  
  
  #---------------- Evaluation objectif
  XB <- data$X %*% B
  Mu <- VectorToMatrix(XB, n,p)
  MuMC <- (Mu + M %*% t(C))
  SC2 <- S %*% t(C*C)
  A <- sqrt(MuMC^2 +SC2)
  expmA <- exp(-A)
  D <- 1 / (1 + expmA)
  E <- D*(1-expmA) / (4 * A)
  
  U <- sum(R*(-log(D)+A-D*A +E*A^2))
  U2 <- sum(R*(Y-D + 2*A*E)*MuMC)
  U3 <-  sum(R*E*(SC2 + MuMC^2))
  U4 <- - 1/2 *sum(M^2 +S - log(S)) + n*q/2
  objective1 <- -U +U2 - U3 + U4   
  
  #T1 <- sum( R* (Y - D) * (MuMC))
  #T2 <- - sum( R* E * ( SC2 + (MuMC - A)^2)) 
  #G <- log(1 + expmA) + A  - D * A
  #T3 <- - sum(R * G) 
  #objective2 <- T1 + T2 + T3 + U4
  
  return(c(-objective1))
}

