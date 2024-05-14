library(Rcpp)

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_BLN_B.cpp"))


#--------------Jeu de données----------------

# Simulation d'un jeu de données simple

n <- 100
p <- 10
d <- 2
q <- 5

X <- matrix(rnorm(n*p*d), nrow = n*p)
X[,1] <- 1
B <- c(-2, 1) #, 10, 2)

XB <- VectorToMatrix(X%*%B, n, p)
Z <- XB  

Prob <- plogis(Z)

Y <- matrix(rbinom(n*p, prob = Prob, size = 1), nrow = n)

vecY <- MatrixToVector(Y)
fit <- glm(vecY ~ -1 + X, family = "binomial")

B.logit <- as.matrix(fit$coefficients)
B.init <- matrix(1,nrow= d, ncol =1)

params <- list(B = B.init)

data <- list(Y = Y,
             R = matrix(1, n, p),
             X = X)

config <- PLNPCA_param()$config_optim



out <- nlopt_optimize_logit(data, params, config)







