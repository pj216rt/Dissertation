library(tidyverse)
library(splines) #for working with b splines
library(mvtnorm) #multivariate normal
library(spam)
library(MASS)
library(MCMCpack) #for the inverse wishart distribution
library(future)
library(future.apply)
library(Matrix)
library(cmdstanr)
library(progressr)

handlers(global = TRUE)

#generate smooth beta functions for the trivariate functional data
beta.gen.fun.tri <- function(grid.size, knots, parameters, deg){
  x <- seq(0, grid.size, length.out = grid.size)
  
  #b spline basis matrix
  bs_basis <- bs(x, degree = deg)  
  
  beta_list <- vector("list", 3)
  
  #going to loop over each dimension
  for (dim in 1:3) {
    random_coeffs_matrix <- matrix(rnorm(ncol(bs_basis)*parameters), 
                                   nrow = ncol(bs_basis), ncol = parameters)
    
    #compute the beta matrix for dimension dim
    beta_matrix <- bs_basis %*% random_coeffs_matrix
    
    #store this in the list
    beta_list[[dim]] <- t(beta_matrix)
  }
  
  #combine all three matrices column wise
  coefs <- do.call(cbind, beta_list)
  
  return(coefs)
}

#function to generate a random covariance matrix
gen.covariance.mat <- function(D, type = c("random", "ar1", "cs", "iid"),
                               rho = 0.6, sigma=1){
  #bunch of cases
  #AR1
  if (type == "ar1") {
    D <- as.matrix(dist(1:D, diag = TRUE, upper = TRUE))
    Sigma <- (rho^D) * sigma^2
  }
  #compound symmetry
  else if (type == "cs") {
    Sigma <- matrix(rho, D, D)
    diag(Sigma) <- 1
    Sigma <- Sigma * sigma^2
  }
  #random, no structure, just PD covariance matrix
  else if (type == "random") {
    A <- matrix(rnorm(D^2), D, D)
    Sigma <- crossprod(A)
  }
  #diagonal matrix
  else if (type == "iid"){
    Sigma <- diag(sigma^2, D)
  }
  return(Sigma)
}

#generate 3 dimensional data, assuming a balanced simulation design
gen.3d.data <- function(num.subj = 50, num.visits = 5,
                        beta.vals, var.x = 1, var.z=1,
                        cov.type = "ar1", rho = 0.5, sigma = 1){
  #accepts num.subj = the number of subjects
  #num.visits = the number of repeated visits/repeated measurements
  #beta.vals = the matrix of true coefficient curves
  #var.x = the variance of the covariates
  #var.z = the variance of the random effects
  
  #error message
  #make sure you supply beta vals
  if (missing(beta.vals)) {
    message("Beta values not supplied.")
  }
  
  #get the grid size and number of parameters from the beta coefficient matrix
  grid.size.3d <- ncol(beta.vals)
  n.params <- nrow(beta.vals)
  
  #compute total observations
  tot.obs <- num.subj*num.visits
  
  #generate the X design matrix
  #each subject receives their own set of covariates
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  
  X$subj = as.factor(X$subj)
  
  #X1 replicates each subject's covariates num.visits times to put it in long format
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  
  #create design matrix dropping intercept and the column subj
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects by multiplying design matrix by beta values
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  
  #create random effect design matrix using subject IDs
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  
  #creates a matrix where each row is a subject, and each column corresponds to
  #a point on the grid.  Basically each row is a set of different random deviations
  #over the grid
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size.3d), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size.3d)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals, adding Gaussian noise for each observation
  #eps <- matrix(rnorm(n = tot.obs * grid.size.3d), nrow = tot.obs, ncol = grid.size.3d)
  Sigma_eps <- gen.covariance.mat(D = grid.size.3d, type = cov.type, rho = rho, sigma = sigma)
  eps <- t(apply(matrix(1:tot.obs), 1, function(i)
    MASS::mvrnorm(1, mu = rep(0, grid.size.3d), Sigma = Sigma_eps)
  ))
  
  #add random effects and epsilon errors to the observed data
  #Yij.true = a true latent trajectory??
  Yij.true <- fixef + ranef
  
  #Yij.obs = add level 1 error
  Yij.obs <- fixef + ranef + eps
  
  #return "observed" data and other stuff in a named list
  output.list <- list()
  output.list$obs <- Yij.obs
  output.list$true <- Yij.true
  output.list$x_design <- X.des
  output.list$z_design <- Z.des
  output.list$raw.data <- X1
  output.list$Sigma_eps <- Sigma_eps
  output.list$cov.type <- cov.type
  output.list$rho <- rho
  output.list$sigma <- sigma
  return(output.list)
  
}

#Goldsmith's modified Gibbs sampler for 3 dimensions
goldsmiths.three.d.gibbs <- function(Y, fixed_design.mat, random_design.mat, 
                                     Kt, N.iter = 1000, N.burn = 200, alpha = .1){
  ##This Gibbs sampler accepts a Y matrix composed of concatenated observations
  ##for the 3 dimensions, a fixed effect design matrix, a random effect design
  ##matrix, the number of knots to be fit for each dimension, as well as overrideable 
  ##defaults for the number of iterations and number of burn in iterations.  One can
  ##also change the alpha term.
  
  
  set.seed(1)
  
  ## fixed and random effect design matrices
  # W.des = model.matrix( fixef.form, data = data)
  # W.des <- W.des[rep(seq_len(nrow(W.des)), length(id)), ]
  # Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  # W.des = as.spam(W.des)
  # print(dim(W.des))
  # Z.des = as.spam(Z.des)
  ## fixed and random effect design matrices
  #W.des = model.matrix( fixef.form, data = data)
  #Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  W.des = as.spam(fixed_design.mat)
  Z.des = as.spam(random_design.mat)
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]/3 #divide by three here
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  Gamma = kronecker(diag(1,3), Theta)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  # not doing anything to the penalty matrix here
  #need to wait until we declare the variance terms first
  
  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))
  
  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  
  Wi = W.des[firstobs, ]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(kronecker(diag(1, I, I), P.mat), diag(1,3)) #kronecker this by a 3*3 matrix
  WIk = kronecker(Wi, diag(1, 3*Kt, 3*Kt)) #double the knots
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP
  
  # initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(Z.des)%*% Z.des, t(Gamma) %*% Gamma)) %*% t(kronecker(Z.des, Gamma)) %*% Y.vec
  bz = matrix(vec.bz, nrow = 3*Kt, ncol = I) #double the knots
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, nrow = 3*Kt, ncol = p) #need to double the knots
  
  Yhat = as.matrix(Z.des %*% t(bz) %*% t(Gamma))
  varhat = var(as.vector(Y - Yhat))
  
  Psi = diag(varhat*IJ, 3*D, 3*D)
  v = IJ
  inv.sig = solve(Psi/v)
  
  Az = I*Kt / 2  #set Az equal to Az1 equal to Az2
  Az1 = I*Kt / 2 
  Az2 = I*Kt / 2
  Bz = sum(diag((t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ])) %*% P.mat %*% t(t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ]))))
  Bz1 = sum(diag((t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ]))))
  Bz2 = sum(diag((t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ]))))
  
  Aw = Kt / 2  #set Aw equal to Aw1
  Aw1 = Kt / 2
  Aw2 = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[1:Kt,u]) %*% P.mat %*% (bw[1:Kt,u])))))
  Bw1 = sapply(1:p, function(u) max(1, sum(diag( t(bw[(Kt+1):(2*Kt),u]) %*% P.mat %*% (bw[(Kt+1):(2*Kt),u])))))
  Bw2 = sapply(1:p, function(u) max(1, sum(diag( t(bw[((2*Kt)+1):(3*Kt),u]) %*% P.mat %*% (bw[((2*Kt)+1):(3*Kt),u])))))
  
  ## matrices to store within-iteration estimates
  BW = array(NA, c(3*Kt, p, N.iter)) #change t0 2*Kt
  BW[,,1] = bw
  BZ = array(NA, c(3*Kt, I, N.iter)) #change to 2*Kt
  BZ[,,1] = bz
  INV.SIG = array(NA, c(3*D, 3*D, N.iter)) #change to 2D
  INV.SIG[,,1] = inv.sig
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW[1,] = lambda.bw = Aw/Bw
  LAMBDA.BW.1 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.1[1,] = lambda.bw.1 = Aw1/Bw1
  LAMBDA.BW.2 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.2[1,] = lambda.bw.2 = Aw2/Bw2
  LAMBDA.BZ = rep(NA, N.iter)
  LAMBDA.BZ[1] = lambda.ranef = Az/Bz
  LAMBDA.BZ.1 = rep(NA, N.iter)
  LAMBDA.BZ.1[1] = lambda.ranef.1 = Az1/Bz1
  LAMBDA.BZ.2 = rep(NA, N.iter)
  LAMBDA.BZ.2[1] = lambda.ranef.2 = Az2/Bz2
  
  y.post = array(NA, dim = c(IJ, 3*D, (N.iter - N.burn))) #change dimensions to 2*D
  
  cat("Beginning Sampler \n")
  pb <- txtProgressBar(min = 0, max = N.iter, initial = 0, style = 3)
  for(i in 1:N.iter){
    setTxtProgressBar(pb,i)
    # if(i %% 100 == 0){ #adding a print statement to tell us where we are
    #   print(i)
    # }
    
    #stick these outside of the subj loop?
    combined <- diag(c(lambda.ranef, lambda.ranef.1, lambda.ranef.2))
    two.P = kronecker(P.mat, combined)  #constructing P \otimes two variances
    
    #sigma_w_k.pre.kron <- c(LAMBDA.BW[i, ], LAMBDA.BW.1[i, ])
    sigma_w_k.pre.kron <- c(lambda.bw, lambda.bw.1, lambda.bw.2)
    #print(sigma_w_k.pre.kron)
    sigma_w_k.pre.kron.1 <- diag(sigma_w_k.pre.kron, nrow = length(sigma_w_k.pre.kron))
    sigma_w_k.post.kron <- kronecker(sigma_w_k.pre.kron.1, P.mat)
    #print(dim(sigma_w_k.post.kron))
    
    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################
    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Gamma))   #change to Gamma
      
      
      #print(dim(two.P))
      sigma = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) +
                      two.P)
      mu = sigma %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) +
                        (two.P) %*% bw %*% t(Wi[subj,]))
      
      bz[,subj] = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = 1)
    }
    ranef.cur = Z.des %*% t(bz) %*% t(Gamma)
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    #error in evaluating the argument 'a' in selecting a method for function 
    #'solve': NAs in argument 5 and 'NAOK = FALSE' (dotCall64)
    
    #maybe here?
    #IIP.1 <- kronecker(kronecker(diag(1, I, I), P.mat), combined)
    IIP.1 <- kronecker(diag(1, I, I), two.P)
    tWIW.1 <- t(WIk) %*% IIP.1 %*% WIk
    
    tWI.1 = t(WIk) %*% IIP.1
    
    sigma = solve(tWIW.1 + sigma_w_k.post.kron)
    #print(determinant(sigma))
    mu = sigma %*% (tWI.1 %*% as.vector(bz))
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = p)
    # 
    beta.cur = t(bw) %*% t(Gamma)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    resid.cur = Y - ranef.cur
    inv.sig = solve(riwish(v + IJ, Psi + t(resid.cur) %*% resid.cur))
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      a.post.1 = Aw1 + Kt/2
      a.post.2 = Aw2 + Kt/2
      b.post = Bw[term] + 1/2 * bw[(1:Kt),term] %*% P.mat %*% bw[(1:Kt),term]
      b.post.1 = Bw1[term] + 1/2 * bw[(Kt+1):(2*Kt),term] %*% P.mat %*% bw[(Kt+1):(2*Kt),term]
      b.post.2 = Bw1[term] + 1/2 * bw[(2*Kt+1):(3*Kt),term] %*% P.mat %*% bw[(2*Kt+1):(3*Kt),term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
      lambda.bw.1[term] = rgamma(1, a.post.1, b.post.1)
      lambda.bw.2[term] = rgamma(1, a.post.2, b.post.2)
    }
    
    ## lambda for random effects
    a.post = Az + I*Kt/2
    a.post.1 = Az1 + I*Kt/2
    a.post.2 = Az2 + I*Kt/2
    b.post = Bz + .5 * sum(sapply(1:I, function(u) (t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) %*% P.mat %*% t(t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) ))
    b.post.1 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) ))
    b.post.2 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) ))
    lambda.ranef = rgamma(1, a.post, b.post)
    lambda.ranef.1 = rgamma(1, a.post.1, b.post.1)
    lambda.ranef.2 = rgamma(1, a.post.2, b.post.2)
    
    ###############################################################
    ## save this iteration's parameters
    ###############################################################
    
    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BW.1[i,] = lambda.bw.1
    LAMBDA.BW.2[i,] = lambda.bw.2
    LAMBDA.BZ[i] = lambda.ranef
    LAMBDA.BZ.1[i] = lambda.ranef.1
    LAMBDA.BZ.2[i] = lambda.ranef.2
    
    if(i > N.burn){
      y.post[,,i - N.burn] = ranef.cur
    }
  }
  close(pb)
  
  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################
  
  #Save the relevant posterior summaries
  
  ## main effects
  beta.post = array(NA, dim = c(p, 3*D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Gamma)
  }
  
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))
  
  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = 3*D)
  for(i in 1:I){
    b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = 3*D)
    for(n in 1:(N.iter - N.burn)){
      b.post[n,] = BZ[,i, n + N.burn] %*% t(Gamma)
    }
    b.pm[i,] = apply(b.post, 2, mean)
  }
  
  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  ranef.pm = Z.des %*% b.pm
  Yhat = apply(y.post, c(1,2), mean)
  
  #returning what we want to return
  ret = list(
    beta.pm = beta.pm,
    beta.LB = beta.LB,
    beta.UB = beta.UB,
    ranef.pm = ranef.pm,
    sig.pm = sig.pm,
    Yhat = Yhat,
    B_W = BW,
    B_Z = BZ,
    INV.SIG = INV.SIG,
    LAMBDA.BW = LAMBDA.BW,
    LAMBDA.BW.1 = LAMBDA.BW.1,
    LAMBDA.BW.2 = LAMBDA.BW.2,
    LAMBDA.BZ = LAMBDA.BZ,
    LAMBDA.BZ.1 = LAMBDA.BZ.1,
    LAMBDA.BZ.2 = LAMBDA.BZ.2,
    y.post = y.post
  )
  
  return(ret)
}

#testing
betas <- beta.gen.fun.tri(grid.size = 25, knots = 20, parameters = 5, deg = 3)
data <- gen.3d.data(beta.vals = betas, cov.type = "random")
#three.d.gibbs.output <- goldsmiths.three.d.gibbs(Y = data$obs, fixed_design.mat = data$x_design,
                                                 #random_design.mat = data$z_design, Kt = 10)

init_fun <- function() {
  list(
    L_Omega_y = diag(1, D),          
    sigma_y   = rep(0.5, D),
    sigma_Z   = rep(0.5, 3),
    sigma_W   = array(0.5, dim = c(p, 3)),
    B_W_raw   = matrix(0, 3 * Kt, p),
    B_Z_raw   = matrix(0, 3 * Kt, I)
  )
}

#compiling test code
test_model <- cmdstan_model("copy_of_gibbs.stan", force_recompile = TRUE)

#set up data
I  <- ncol(data$z_design)
IJ <- nrow(data$obs)
p <- ncol(data$x_design)
D  <- ncol(data$obs)
D0 <- D/3
Kt <- 10
J <- 5
alpha <- 0.01

Theta <- splines::bs(1:D0, df = Kt, intercept = TRUE, degree = 3)
THETA <- kronecker(diag(3), Theta)

#penalty matrix
diff2  <- matrix(rep(c(1,-2,1, rep(0, D0-2)), D0-2)[1:((D0-2)*D0)], D0-2, D0, byrow=TRUE)
P0     <- crossprod(Theta)
P2     <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
P.mat  <- alpha*P0 + (1 - alpha)*P2
PenMatInv <- solve(P.mat)
PenMatInv_reg <- PenMatInv + 1e-6 * diag(Kt)
scale_pen <- mean(diag(PenMatInv_reg))
PenMatInv_scaled <- PenMatInv_reg / scale_pen

SUBJ = factor(apply(data$z_design %*% 1:I, 1, sum))

firstobs <- sapply(unique(SUBJ), function(s) which(SUBJ==s)[1])
Wi <- data$x_design[firstobs, ]
X <- Wi
Z <- data$z_design

Y <- data$obs

#compile stan data
stan_data <- list(
  I  = I,
  J  = J,
  IJ = IJ,
  D  = D,
  p  = p,
  Kt = Kt,
  
  Y = Y,          
  X = X,          
  subj_id = max.col(data$z_design, ties.method = "first"),
  
  THETA = THETA,      
  PenMatInv = PenMatInv_scaled
  #PenMatInv = diag(Kt)
)

fit_stan_test <- test_model$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 200,
    iter_sampling = 100,
    seed = 123,
    refresh = 5,
    init = init_fun
  )


#compiling  old stan code
old_stan <- cmdstan_model("corr_rand_effects_model.stan")


# #testing new stan code
# new_stan <- cmdstan_model("indw_z_corr_Sigma.stan")

# I  <- ncol(data$z_design)                 #subjects
# D  <- ncol(data$obs) / 3                  #per-dimension grid
# D3 <- 3L * D
# IJ <- nrow(data$obs)
# p  <- ncol(data$x_design)
# Kt <- 10
# alpha <- 0.01
# 
# #basis
# Theta <- splines::bs(1:D, df = Kt, intercept = TRUE, degree = 3)  # D x Kt
# 
# #penalty matrix
# diff2  <- matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow=TRUE)
# P0     <- crossprod(Theta)
# P2     <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
# P.mat  <- alpha * P0 + (1 - alpha) * P2
# 
# SUBJ = factor(apply(data$z_design %*% 1:dim(data$z_design)[2], 1, sum))
# 
# #first observation
# firstobs = rep(NA, length(unique(SUBJ)))
# for(i in 1:length(unique(SUBJ))){
#   firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
# }
# 
# Wi = data$x_design[firstobs, ]
# 
# #Gamma
# THETA3 <- kronecker(diag(3), Theta)                            # (D3 x 3Kt)
# 
# a_Z <- rep(2, 3); b_Z <- rep(1, 3)
# a_W <- matrix(2, nrow = 3, ncol = p)
# b_W <- matrix(1, nrow = 3, ncol = p)
# 
# #subject index
# subj_id <- max.col(data$z_design, ties.method = "first")
# 
# #prepare data for stan
# stan_data <- list(
#   I = I, J = 5,
#   IJ = IJ,
#   D = D,                  # per-dimension grid (keep, if your Stan uses it)
#   D3 = D3,                # total stacked grid (3D) -> NEW
#   p = p,
#   Kt = Kt,
#   Y = as.matrix(data$obs),             # IJ x D3
#   #X = as.matrix(data$x_design),        # I x p
#   X = Wi,
#   #Z = as.matrix(data$z_design),        # IJ x I
#   THETA = THETA3,                     # D3 x 3Kt (Gamma)
#   PenMatInv = solve(P.mat),                      # preferred with multi_normal_prec
#   # PenMatInv = symmetrized_solve_if_your_model_needs_it,
#   
#   corr_code = 3,
#   a_Z = a_Z, b_Z = b_Z,
#   a_W = a_W, b_W = b_W,
#   
#   subj_id = as.integer(subj_id)
# )
# 
# fit_stan_test <- new_stan$sample(
#   data = stan_data,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 2000,
#   iter_sampling = 1000,
#   seed = 123,
#   refresh = 5,
#   adapt_delta = 0.98,
#   max_treedepth = 12
# )


cov_codes <- c(iid = 1, cs = 2, ar1 = 3)
#stan_model_structured <- cmdstan_model("fos_correlation_structure_random_effects.stan")

#wrapper function to fit and time the Gibbs sampler
fit_gibbs_wrapper <- function(data, Kt = 10, iter = 500, burn = 100) {
  t0 <- Sys.time()
  fit <- goldsmiths.three.d.gibbs(
    Y = data$obs,
    fixed_design.mat = data$x_design,
    random_design.mat = data$z_design,
    Kt = Kt, N.iter = iter, N.burn = burn
  )
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  list(fit = fit, time_sec = sec)
}

#wrapper function to fit and time the STAN sampler
fit_stan_wrapper <- function(data, Kt, iter_warmup = 1500, iter_sampling = 1000,
                             seed = 1234, adapt_delta = 0.98, max_treedepth = 10){
  
  #compiling stan model
  #mod <- cmdstan_model("corr_rand_effects_model.stan")
  mod <- cmdstanr::cmdstan_model(exe_file = model_exe_path)
  
  #set up data
  I  <- ncol(data$z_design)
  IJ <- nrow(data$obs)
  p <- ncol(data$x_design)
  D  <- ncol(data$obs)
  D0 <- D/3
  #Kt <- 10
  J <- 5
  alpha <- 0.01
  
  Theta <- splines::bs(1:D0, df = Kt, intercept = TRUE, degree = 3)
  THETA <- kronecker(diag(3), Theta)
  
  #penalty matrix
  diff2  <- matrix(rep(c(1,-2,1, rep(0, D0-2)), D0-2)[1:((D0-2)*D0)], D0-2, D0, byrow=TRUE)
  P0     <- crossprod(Theta)
  P2     <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat  <- alpha*P0 + (1 - alpha)*P2
  PenMatInv <- solve(P.mat)
  
  SUBJ = factor(apply(data$z_design %*% 1:I, 1, sum))
  
  firstobs <- sapply(unique(SUBJ), function(s) which(SUBJ==s)[1])
  Wi <- data$x_design[firstobs, ]
  X <- Wi
  Z <- data$z_design
  
  Y <- data$obs
  
  #compile stan data
  stan_data <- list(
    I  = I,
    J  = J,
    IJ = IJ,
    D  = D,
    p  = p,
    Kt = Kt,
    
    Y = Y,          
    X = X,          
    subj_id = max.col(data$z_design, ties.method = "first"),
    
    THETA = THETA,      
    PenMatInv = PenMatInv
  )
  
  #run stan
  t0 <- Sys.time()
  fit <- mod$sample(
    data = stan_data,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = 100,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  #need to extract the posterior draws
  BETA_draws <- fit$draws("BETA", format = "draws_matrix")
  
  S <- nrow(BETA_draws)
  K <- ncol(BETA_draws)
  
  BETA_pm_vec <- colMeans(BETA_draws)
  BETA_pm <- matrix(BETA_pm_vec, nrow = p, ncol = D, byrow = TRUE)
  BETA_draws <- array(BETA_draws, dim = c(S, p, D))
  #BETA_pm <- apply(BETA_draws, c(2,3), mean)
  
  lev1_sigma_draws <- fit$draws("lev1_sigma")  
  Sigma_pm <- mean(lev1_sigma_draws)^2*diag(D)
  
  return(list(
    beta_pm     = BETA_pm,        
    BETA_draws  = BETA_draws,     
    Sigma_pm    = Sigma_pm,       
    time_sec    = sec,
    stan_fit    = fit             
  ))
}

#need to design simulation study
sim_study_compare <- function(
    R = 1,
    subj.levels = c(25),
    true.cov.types = c("iid","cs","ar1","random"),
    fit.cov.types  = c("iid","cs","ar1"),
    rho.levels = c(0.5),
    sigma_z_levels = c(1),
    sigma1_levels = c(1),
    n_cores = max(1, parallel::detectCores() - 1),
    iter = 600, burn = 100,
    K.levels = c(10)
){
  
  #precent cmdstanr from recompiling
  Sys.setenv(CMDSTANR_NO_AUTO_RECOMPILE = "true")
  
  mod <- cmdstanr::cmdstan_model("corr_rand_effects_model.stan", compile = TRUE)
  mod_path <- mod$exe_file()
  
  library(progressr)
  handlers("progress")
  
  set.seed(2025)
  plan(multisession, workers = n_cores)
  
  #commenting out the fitted covarianace, not doing anything yet.  
  design <- expand.grid(rep = 1:R,
                        subj = subj.levels,
                        true.cov = true.cov.types,
                        #fit.cov = fit.cov.types,
                        rho = rho.levels,
                        sigmaz = sigma_z_levels,
                        sigma1 = sigma1_levels,
                        K = K.levels,
                        KEEP.OUT.ATTRS = FALSE)
  
  pr <- progressor(along = 1:nrow(design))
  
  #constant beta for all runs
  beta_true <- beta.gen.fun.tri(grid.size = 45, knots = 10, parameters = 5, deg = 3)
  
  run_one <- function(ix){
    library(MASS)
    library(spam)
    library(coda)
    library(tibble)
    library(dplyr)
    library(splines)
    
    logfile <- "test_log.txt"
    write(sprintf("[%s] Starting run_one(%d)\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    #implement progress
    pr()
    #defining the number of parameters.  Future apply creates new R sessions
    #need everything needed in run_one to be define within run_one.  
    p <- nrow(beta_true)
    D <- ncol(beta_true)
    
    #getting run specific parameters
    nsub <- design$subj[ix]
    true_cov <- design$true.cov[ix]
    #fit_cov  <- design$fit.cov[ix]
    rho_val  <- design$rho[ix]
    sigz_val <- design$sigma1[ix]
    sig1_val <- design$sigma1[ix]
    K_val <- design$K[ix]
    
    write(sprintf("[%s] Generating data for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    #simulate one dataset
    data <- gen.3d.data(num.subj = nsub, num.visits = 5,
                        beta.vals = beta_true, var.z = sigz_val,
                        cov.type = true_cov, rho = rho_val, sigma = sig1_val)
    
    write(sprintf("[%s] Finished data generation for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    #Gibbs simulation study
    write(sprintf("[%s] Starting Gibbs for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    g <- fit_gibbs_wrapper(data, Kt = K_val, iter = iter, burn = burn)
    
    write(sprintf("[%s] Finished Gibbs for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    #Metrics for Gibbs
    beta_hat_g <- g$fit$beta.pm
    IMSE_g <- mean((beta_hat_g - beta_true)^2)
    Coverage_g <- mean((beta_true >= g$fit$beta.LB) & (beta_true <= g$fit$beta.UB))
    SigmaErr_g <- norm(g$fit$sig.pm - data$Sigma_eps, "F") / norm(data$Sigma_eps, "F")
    
    #ess
    ess_beta_g <- coda::effectiveSize(g$fit$B_W[1,1,])
    
    #stan
    write(sprintf("[%s] Starting Stan for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    s <- fit_stan_wrapper(data, Kt = K_val, iter_warmup = burn,
                          iter_sampling = iter-burn)
    
    write(sprintf("[%s] Finished Stan for row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)
    
    #getting beta hat
    beta_hat_s <- s$beta_pm
    beta_hat_s <- matrix(beta_hat_s, nrow = p, ncol = D, byrow = TRUE)
    
    BETA_draws <- s$BETA_draws
    IMSE_s <- mean((beta_hat_s - beta_true)^2)
    BETA_LB <- apply(BETA_draws, c(2,3), quantile, 0.025)
    BETA_UB <- apply(BETA_draws, c(2,3), quantile, 0.975)
    Coverage_s <- mean((beta_true >= BETA_LB) & (beta_true <= BETA_UB))
    SigmaErr_s <- norm(s$Sigma_pm - data$Sigma_eps, "F") / norm(data$Sigma_eps, "F")
    
    #with just 1 chain, R hat will be useless.  ESS is still valid though
    summ <- s$stan_fit$summary(variables = "BETA[1,1]")
    ess_beta_s <- summ$ess_bulk
    
    write(sprintf("[%s] Finished row %d\n", Sys.time(), ix),
          file = logfile, append = TRUE)

    #compiling the each row
    out <- tibble::tibble(
      rep = design$rep[ix],
      subj = nsub,
      Knots = K_val,
      true.cov = true_cov,
      sig.z = sigz_val,
      sig.lev1 = sig1_val,
      #fit.cov  = fit_cov,
      rho = rho_val,
      
      #gibbs
      IMSE_gibbs     = IMSE_g,
      Coverage_gibbs = Coverage_g,
      SigmaErr_gibbs = SigmaErr_g,
      ESS_gibbs = ess_beta_g,
      Time_gibbs = g$time_sec,
      
      #stan
      IMSE_stan = IMSE_s,
      Coverage_stan = Coverage_s,
      SigmaErr_stan = SigmaErr_s,
      ESS_stan  = ess_beta_s,
      Time_stan  = s$time_sec
    )
    
    #progress bar
    #write(sprintf("Finished row %d\n", ix), file = "progress.log", append = TRUE)
    
    return(out)
  }
  
  model_exe_path <- mod_path 
  
  #run this in parallel across multiple cores
  res_list <- future_lapply(1:nrow(design), run_one, future.seed = TRUE, future.stdout = TRUE,
                            future.globals = TRUE)
  dplyr::bind_rows(res_list)
}


results <- sim_study_compare(R=5, n_cores = 6, iter = 1200, burn = 200, K.levels = c(10, 20), sigma1_levels = c(0.5,1),
                             sigma_z_levels = c(0.1, 0.5))

#print(results)

#save csv 
write.csv(results, "results_sim.csv")
