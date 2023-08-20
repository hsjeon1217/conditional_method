#############################################################################################################################
# library(sigmoid)
# https://rdrr.io/cran/pROC/man/auc.html # library(pROC)

# Motivated from R_functions_0929.R
#############################################################################################################################

pkgs = c("pROC","rpart","Inference","Rcpp","parallel","doParallel","foreach","pipeR","tidyr","dplyr","stats","invgamma","ggplot2","latex2exp","FDRreg")
all(sapply(pkgs , require, character.only = T))
source("Q_value_sources.R")

# Step 0. Data Generation Y & P-values:
# changed from Data.0928

# Data parameters : m, n, Edelta
# The function generating x : fct_pi0

# Updated from Data.0127
Data.0226 = function( m, n, Edelta,
                      fct_pi0 ){

  # log10 gene length x :
  x = rnorm(m, 4, 0.5);
  # x = runif(m, 3,5)

  # Allocate which one is the true DEGs :
  pi0.x = sapply(x, function(.){ fct_pi0(.) })
  true.degs = which( sapply(pi0.x, function(pi0){ as.logical(rbinom(1, size = 1, prob = (1-pi0))) }))

  # Data Generation # intercept = 0 <- Actually not important.
  v_sd = rinvchisq(n = m, df = 5, ncp = 0); # mean = 1/3
  Y.ctrl  <- (matrix( rnorm((n*m), 0, 1), nrow = m, ncol = n ) * v_sd) # m x n
  Y.trt   <- (matrix( rnorm((n*m), 0, 1), nrow = m, ncol = n ) * v_sd) # m x n

  # For the first m1 genes, add fixed effect :
  m1 = length(true.degs); deltas = rnorm(m1, Edelta, sd = 0.02);
  Y.trt[true.degs,] = Y.trt[true.degs,] + deltas;

  # The order of Y and X is trt -> ctrl.
  Y = cbind(Y.trt, Y.ctrl);
  X = model.matrix(~ 0 + as.factor(c(rep(1,n), rep(2,n))));
  colnames(X) = c("trt","ctrl");

  # Traditional T-test from MSE.generate.0521
  df = 2*(n - 1);
  b = solve(t(X) %*% X) %*% t(X) %*% t(Y);
  Yhat = t(X %*% b);
  MSE = apply((Y - Yhat)^2,1,sum)/df;
  mu = t(b); colnames(mu) = colnames(X);

  p_df = data.frame( MSE = MSE, mu) %>%
    dplyr::mutate( Tstat = (trt - ctrl)/sqrt(2*MSE/n),
            p = 2*pt(abs(Tstat), df = df, lower.tail = F))

  # Code check using t.test package:
  # p_vals = sapply(1:m, function(i){ p = t.test(Y[i,1:n], Y[i,-c(1:n)])$p.value})
  # plot(p_vals, p_df$p)

  t.score = p_df$Tstat
  
  return(list( true.degs = true.degs,
               data = data.frame( orig.id = 1:m,
                                  x = x, p = p_df$p ),
               real.pi0.x = pi0.x,
               t.score = t.score
        ))

}



#############################################################################################################################

# Step2. Based on the group info in step1, 


Proposed.Q_0929 = function(data){

  # Input : data.frame(orig.id, f_x, p, pi0) # f_x : cond.pi0
  m = nrow(data)
  
  # Ordered by p.tilde VVV
  temp = data %>% 
    dplyr::mutate(p.tilde = f_x*p) %>% 
    dplyr::arrange(p.tilde) 
  
  v.pi0 = dplyr::pull(temp, "pi0");
  v.f_x = dplyr::pull(temp, "f_x");
  v.p = dplyr::pull(temp, "p"); 
  v.p.tilde = dplyr::pull(temp, "p.tilde");
  
  cst = sum(v.pi0/v.f_x) # sum(pi0/cond.pi0)
  prop.fdr = c(v.p.tilde*cst)/c(1:m) # R(t) = 1:m
  
  # check the q-value VVVVV
  prop.q = rev(cummin(rev(prop.fdr)))
  
  # Ordered by orig.id again. VVV
  temp %>% 
    cbind.data.frame(prop.q = prop.q) %>%
    dplyr::arrange(orig.id) %>% 
    dplyr::select(prop.q) %>%
    return(.)
  
}

#############################################################################################################################

# Wrong code??
# # Motivated from Local.pi0.estimate.1005 - Almost the same. 
# Local.pi0.estimate.0127 = function(x.vec, p.vec, neigh.m = 100, cond.alpha = 0.05, B = 20){
#   
#   # Settings # dist.mat0 = as.matrix(dist(x.vec)) # all(dist.mat0 == dist.mat)
#   m = length(p.vec);
#   dist.mat = abs(matrix(x.vec, nrow = m, ncol = m, byrow = T) - x.vec)  
#   
#   neigh.ids.df <- sapply(1:m, function(i){
#     ids = sort(order(dist.mat[i,])[1:neigh.m])
#   }) %>% data.frame(.)
#   rm(dist.mat); gc();  
#   
#   neigh.p.df = sapply(neigh.ids.df, function(ids){p.vec[ids]}) %>% data.frame(.)
#   
#   # p = cbind(p1, p2)
#   k = ncol(neigh.p.df)
#   bin = seq(0, 1, 1/B)
#   
#   # bin.counts and tail.means
#   bin.info <- apply(t(sapply(bin, function(b){
#     counts = apply(neigh.p.df,2,function(.){sum(.<=b)})
#     return(c(counts))
#   })), 2, diff)
#   
#   # bin.info <- sapply(neigh.p.df, function(p){ unlist(table(cut(p, bin))) })
#   # all(bin.info == bin.info2)
#   tail.info = apply(bin.info, 2, function(.){rev(cumsum(rev(.))/(1:B))}) 
#   
#   indices = sapply(1:k, function(i){
#     min(which(bin.info[,i]<=tail.info[,i]))
#   })
#   
#   cutoffs = bin[-(B+1)][indices] # V treating zero...?
#   
#   # estimate.m0(p.vec[neigh.ids.df[,2]]) == sum(neigh.p.df[,2] > cutoffs[2])/(1-cutoffs[2])
#   
#   pi0.x = sapply(1:m, function(i){
#     sum(neigh.p.df[,i] >= cutoffs[i])/(1-cutoffs[i])/neigh.m
#   })
#   
#   cond.pi0.x = sapply(1:m, function(i){
#     sum(pi0.x[neigh.ids.df[,i]])*cond.alpha/sum(neigh.p.df[,i] <= cond.alpha)
#   })
#   
#   rm(neigh.p.df, neigh.ids.df); gc();
#   
#   return(data.frame( est.pi0.x = pi0.x, 
#                      cond.pi0.x = cond.pi0.x))
#   
# }

# Motivated from Local.pi0.estimate.1018 # Almost the same.
Local.pi0.estimate.0127 = function(x.vec, p.vec, neigh.m = 100, cond.alpha = 0.05, B = 20){
  
  ### Settings # dist.mat0 = as.matrix(dist(x.vec)) # all(dist.mat0 == dist.mat)
  m = length(p.vec);
  
  dist.mat = abs(matrix(x.vec, nrow = m, ncol = m, byrow = T) - x.vec)
  
  # cf. sort.neigh.ids.df = apply(neigh.ids.df, 2, sort)
  neigh.ids.df <- sapply(1:m, function(i){
    ids = order(dist.mat[i,])[1:neigh.m]
  }) %>% data.frame(.)
  
  # rm(dist.mat); gc();
  neigh.p.df = sapply(neigh.ids.df, function(ids){p.vec[ids]}) %>% data.frame(.)
  
  ### Cutoff - Calculation :
  k = ncol(neigh.p.df)
  bin = seq(0, 1, 1/B)
  
  # bin.counts and tail.means
  bin.info <- apply(t(sapply(bin, function(b){
    counts = apply(neigh.p.df,2,function(.){sum(.<=b)})
    return(c(counts))
  })), 2, diff)
  
  # bin.info <- sapply(neigh.p.df, function(p){ unlist(table(cut(p, bin))) })
  # all(bin.info == bin.info2)
  tail.info = apply(bin.info, 2, function(.){rev(cumsum(rev(.))/(1:B))}) 
  
  indices = sapply(1:k, function(i){
    min(which(bin.info[,i]<=tail.info[,i]))
  })
  
  cutoffs = bin[-(B+1)][indices] # V treating zero...?
  
  # estimate.m0(p.vec[neigh.ids.df[,2]]) == sum(neigh.p.df[,2] > cutoffs[2])/(1-cutoffs[2])
  
  # 1. Based on neighborhood :
  pi0 = sapply(1:m, function(i){
    sum(neigh.p.df[,i] >= cutoffs[i])/(1-cutoffs[i])/neigh.m
  })
  
  prop.leq.alpha = unlist(apply(neigh.p.df, 2, function(.){sum(.<=cond.alpha)/neigh.m}))
  
  cond.pi0 =  pi0/prop.leq.alpha*cond.alpha
  cond.pi0[cond.pi0>1] = 1 # VVV
  
  # neigh.q = sapply(neigh.p.df, function(.){ q = jabes.q(.)[1];})
  gc();
  
  rm(dist.mat, neigh.ids.df, neigh.p.df); gc();
  
  # return(list( pi0.df = data.frame( Local.pi0.x = pi0, Local.cond.pi0.x = cond.pi0),
  #              neigh.q = neigh.q))
  
  return(data.frame( est.pi0.x = pi0, cond.pi0.x = cond.pi0))
  
}

#############################################################################################################################

# Data.0127 = function( m, n, Edelta, 
#                       fct_pi0 ){
#   
#   # log10 gene length x :
#   x = rnorm(m, 4, 0.5);
#   # x = runif(m, 3,5)
#   
#   # Allocate which one is the true DEGs :
#   pi0.x = sapply(x, function(.){ fct_pi0(.) }) 
#   true.degs = which( sapply(pi0.x, function(pi0){ as.logical(rbinom(1, size = 1, prob = (1-pi0))) }))
#   
#   # Data Generation # intercept = 0 <- Actually not important.
#   v_sd = rinvchisq(n = m, df = 5, ncp = 0); # mean = 1/3  
#   Y.ctrl  <- (matrix( rnorm((n*m), 0, 1), nrow = m, ncol = n ) * v_sd) # m x n
#   Y.trt   <- (matrix( rnorm((n*m), 0, 1), nrow = m, ncol = n ) * v_sd) # m x n
#   
#   # For the first m1 genes, add fixed effect :
#   m1 = length(true.degs); deltas = rnorm(m1, Edelta, sd = 0.02);
#   Y.trt[true.degs,] = Y.trt[true.degs,] + deltas;
#   
#   # The order of Y and X is trt -> ctrl.
#   Y = cbind(Y.trt, Y.ctrl);
#   X = model.matrix(~ 0 + as.factor(c(rep(1,n), rep(2,n))));
#   colnames(X) = c("trt","ctrl");
#   
#   # Traditional T-test from MSE.generate.0521
#   df = 2*(n - 1);
#   b = solve(t(X) %*% X) %*% t(X) %*% t(Y); 
#   Yhat = t(X %*% b);
#   MSE = apply((Y - Yhat)^2,1,sum)/df;
#   mu = t(b); colnames(mu) = colnames(X);
#   
#   p_df = data.frame( MSE = MSE, mu) %>% 
#     dplyr::mutate( Tstat = (trt - ctrl)/sqrt(2*MSE/n), 
#             p = 2*pt(abs(Tstat), df = df, lower.tail = F))
# 
#   # Code check using t.test package:  
#   # p_vals = sapply(1:m, function(i){ p = t.test(Y[i,1:n], Y[i,-c(1:n)])$p.value})
#   # plot(p_vals, p_df$p)
# 
#   return(list( true.degs = true.degs,
#                data = data.frame( orig.id = 1:m, 
#                                   x = x, p = p_df$p ),
#                real.pi0.x = pi0.x
#         ))
#   
# }