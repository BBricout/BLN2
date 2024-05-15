library(Rcpp)
library(RcppArmadillo)
library(nloptr)
library(missForest)

#---------------------------------------------------------
#------------Simulation d'un jeu de données-------------------
#----------------------------------------------------------------------------
n <- 1000
p <- 20 
d <- 2
q <- 1
list_dim <- list(d=d,n=n,p=p,q=q)

X <- matrix(rnorm(n*p*d), nrow = n*p,ncol=d)
X[,1] <- 1
B <- c(-2, 1, 3, 2, 5, -3)[1:d]
rho <- 2

W <- matrix(rnorm(n*q), nrow = n)
C <- matrix(rho*rnorm(p*q), nrow = p)
XB <- VectorToMatrix(X%*%B, n, p)
Z <- XB + W%*%t(C)    
Prob <- plogis(Z)
Y <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)
data <- list(Y = Y, R = matrix(1, n, p), X = X)

Y.na30 <- missForest::prodNA(Y, 0.3)
data.na30 <- list(Y = Y.na30, R = ifelse(is.na(Y.na30), 0, 1), X = X)

#--------------------------------
########  INITIALISATION
#---------------------------------

params.list.true <- list(B=B,C=C,M = matrix(W,n,q),Sd=sqrt(matrix(0.1,n,q)))
params.list.init <- Init_BLN(data)
params.vec.init <- paramListToVector(params.list.init,list_dim)


Obj_init <- objective(params.vec = params.vec.init, data,list_dim= list(d=d,n=n,p=p,q=q))
Obj_true <- objective(params.vec = paramListToVector(params.list.true,list_dim), data,list_dim= list(d=d,n=n,p=p,q=q))

#################  ESTIMATION only B 
res <- optim(par = params.list.init$B, fn=objective,  data = data,list_dim=list_dim ,params.list.true = params.list.true,estim=c('B'))
cbind(params.list.init$B,res$par,params.list.true$B)
c(Obj_init, res$value,Obj_true)


Zhat <- VectorToMatrix(X%*%params.list.estim$B,n,p)+W%*%t(params.list.true$C)
Zhat.init <- VectorToMatrix(X%*%params.list.init$B,n,p)+W%*%t(params.list.true$C)
plot(Z,Zhat,ylab='true',xlab='estim')
points(Z,Zhat.init,col=('red'))

plot(cov(Y),params.list.estim$C%*%t(params.list.estim$C))

############  ESTIMATION only C 
params.list.init$B <- params.list.true$B
params.init.vec <- paramListToVector(params.list.true,list_dim,estim='C')

res <- optim(par = params.init.vec, fn=objective,  data = data,list_dim=list_dim ,params.list.true = params.list.true,estim=c('C'))
params.list.estim <- VectorToparamsList(res$par,list_dim,estim='C',params.list.true = params.list.true)
c(res$value,Obj_true)

CC.estim <- params.list.estim$C%*%t(params.list.estim$C)
CC.true <- params.list.true$C%*%t(params.list.true$C)
plot(CC.estim,CC.true,ylab='true',xlab='estim')
abline(a=0,b=1)
####################################"" 

################# Estimation B and C. (to do) 
XB + W%*%t(C)    

abline(a=0,b=1)


result <- nloptr(x0 = paramsInit.vec, eval_f = objective, 
                 opts = list("algorithm" = "NLOPT_GN_ISRES", "xtol_rel" = 1.0e-8),
                 data = data,list_dim=list(d=d,n=n,p=p,q=q))
result$x0



