
DDE.proportion.0228 = function( x.vec, p.vec, oob.ids, 
                                FDR.level = 0.05, 
                                alpha.vec = seq(0.01, 0.2, by = 0.01),
                                neigh.m = 2000, B = 20 ){
  
  # Step 0 : Input variables :
  m = length(x.vec); inb.m = (m - length(oob.ids));
  temp.x = c(x.vec[-oob.ids], x.vec[oob.ids]);
  temp.p = c(p.vec[-oob.ids], p.vec[oob.ids]);
  
  # m x inb.m distance matrix :
  dist.mat = abs(matrix(temp.x, nrow = m, ncol = m, byrow = T) - temp.x)[,1:inb.m];
  # dim(dist.mat)
  
  # Step 1:
  # cf. sort.neigh.ids.df = apply(neigh.ids.df, 2, sort)
  neigh.ids.df <- sapply(1:m, function(i){
    ids = order(dist.mat[i,])[1:neigh.m]
  }) %>% data.frame(.)
  # dim(neigh.ids.df)
  
  # rm(dist.mat); gc(); VVVVVVVVVVVVVVVVVVVVVVVVVV p.vec -> temp.p here VVVVVV
  neigh.p.df = sapply(neigh.ids.df, function(ids){temp.p[ids]}) %>% data.frame(.)
  # dim(neigh.p.df) 2000 x m
  
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
  
  # 2-1: Based on neighborhood :
  pi0.vec = sapply(1:m, function(i){
    sum(neigh.p.df[,i] >= cutoffs[i])/(1-cutoffs[i])/neigh.m
  })
  
  # 2-2: Prop(pi <= alpha) -> "prop.df"
  colnames(neigh.p.df) = 1:m
  neigh.p.df %>% 
    gather(., id, neigh.p, colnames(neigh.p.df), factor_key = T) %>% 
    group_by(id) %>%
    mutate(cuts = cut(neigh.p, c(0, alpha.vec))) %>% 
    group_by(id, cuts) %>%
    summarise(n = n(), .groups = "drop") %>%
    spread(., cuts, n, drop = F) -> n.df

  # colnames(n.df) == 1:m
  n.df = n.df[,-ncol(n.df)]
  n.df[is.na(n.df)] = 0; n.df = data.frame(n.df); n.df = t(apply(n.df[,-1], 1, cumsum));
  prop.df = n.df/neigh.m;
  colnames(prop.df) = alpha.vec[-length(alpha.vec)]; # be careful about column names VVVVV
  
  # Find t.star for a given alpha.
  # Using Proposed.Q_0929(data) # data.frame(orig.id, f_x, p, pi0) # f_x : cond.pi0
  source.df = data.frame(orig.id = 1:m, p = temp.p, pi0 = pi0.vec, prop.df) %>% 
                    gather(., alpha, prop, paste0("X",colnames(prop.df)), factor_key = T) %>% 
                    mutate( alpha = as.numeric(as.character(gsub("X(.*)","\\1",alpha))),
                            f_x = alpha * pi0/prop )
  # conditional probability
  source.df$f_x[source.df$f_x>1] = 1 
  
  
  inb.source.df = source.df %>% filter(orig.id <= inb.m) %>% dplyr::select(alpha, orig.id, f_x, p, pi0)
  t.star.vec = sapply(sort(unique(inb.source.df$alpha)), function(alpha){
    # print(alpha)
    data = inb.source.df[inb.source.df$alpha == alpha, -1]
    q.values = Proposed.Q_0929(data)
    # Check here !!!!
    temp.df = data[which(q.values$prop.q <= FDR.level),]
    t.star = max(temp.df$p * temp.df$f_x)
    # print(t.star)
    return(t.star)
  })
  
  # Using the t.star.vec, make inference on oob.data
  source.df %>% filter(orig.id > inb.m) %>% 
    left_join(., data.frame(alpha = sort(unique(inb.source.df$alpha)), t.star = t.star.vec), by = "alpha") %>%
    mutate(p.tilde = p*f_x) %>% 
    mutate(dde = p.tilde <= t.star) %>% group_by(alpha) %>% 
    summarise(., sum(dde)/inb.m, .groups = "drop") -> DDE.prop

  return(DDE.prop)
  
}



# nfold = 3
CV.alpha.0228 = function(x.vec, p.vec, FDR.level, neigh.m, nfold, B = 20){

# input variable :
m = length(x.vec);
# candidate alphas :
alpha.vec = c(seq(0.01, 0.2, by = 0.01), 1);

# k - fold cross-validation :
folds = sample(rep(seq(nfold), length = m));

lapply(1:nfold, function(f){
# print(paste0("cross-validation fold = ", f))
oob.ids = sort(which(folds == f))
DDE.prop.vec = DDE.proportion.0228( x.vec, p.vec, oob.ids, FDR.level, alpha.vec, neigh.m, B)

}) %>% Reduce("+",.)/nfold -> DDE.prop.vec

colnames(DDE.prop.vec)[2] = "DDE.prop"
DDE.prop = DDE.prop.vec$DDE.prop

# print(DDE.prop.vec$DDE.prop)
# alpha.star = alpha.vec[which.max(DDE.prop.vec$DDE.prop)]

return(list(alpha = alpha.vec, DDE.prop = DDE.prop))

}


CV.alpha.0324 = function(alpha.vec, x.vec, p.vec, FDR.level, neigh.m, nfold, B = 20){
  
  # input variable :
  m = length(x.vec);
  # candidate alphas :
  # alpha.vec = c(seq(0.01, 0.2, by = 0.01), 1);
  
  # k - fold cross-validation :
  folds = sample(rep(seq(nfold), length = m));
  
  lapply(1:nfold, function(f){
    # print(paste0("cross-validation fold = ", f))
    oob.ids = sort(which(folds == f))
    DDE.prop.vec = DDE.proportion.0228( x.vec, p.vec, oob.ids, FDR.level, alpha.vec, neigh.m, B)
    
  }) %>% Reduce("+",.)/nfold -> DDE.prop.vec
  
  colnames(DDE.prop.vec)[2] = "DDE.prop"
  DDE.prop = DDE.prop.vec$DDE.prop
  
  # print(DDE.prop.vec$DDE.prop)
  # alpha.star = alpha.vec[which.max(DDE.prop.vec$DDE.prop)]
  
  return(list(alpha = alpha.vec, DDE.prop = DDE.prop))
  
}


# x.vec = log10(p_df$GeneLength); p.vec = p_df$pB;
# # Inference based on proposed method :
# # Step 1. Cross-Validation to maximize the DDE genes number :
# FDR.level = 0.2; neigh.m = 2000; nfold = 10;
# cv.alpha = CV.alpha.0228(x.vec, p.vec, FDR.level, neigh.m, nfold, B = 20)
# 
# cv.alpha = cv.alpha$alpha[-length(cv.alpha$alpha)][which.max(cv.alpha$DDE.prop)]
# 
# plot(cv.alpha$alpha[-length(cv.alpha$alpha)], cv.alpha$DDE.prop)
# 
# # Step 2. Based on the alpha star, q-values are generated.
# # data = data.frame(orig.id, f_x, p, pi0) 
# # Local.pi0.estimate.0127 -> return(data.frame( est.pi0.x = pi0, cond.pi0.x = cond.pi0))
# 
# cv.source.df = Local.pi0.estimate.0127(x.vec, p.vec, neigh.m, cond.alpha = cv.alpha, B = 20)
# cv.q.values = Proposed.Q_0929(data = data.frame( orig.id = 1:m, 
#                                                  f_x = cv.source.df$cond.pi0.x, 
#                                                  p = p.vec, 
#                                                  pi0 = cv.source.df$est.pi0.x ))
# 
# # default setting 1 alpha = 0.05
# default1.source.df = Local.pi0.estimate.0127(x.vec, p.vec, neigh.m, cond.alpha = 0.05, B = 20)
# default1.q.values = Proposed.Q_0929(data = data.frame(orig.id = 1:m, 
#                                                       f_x = default1.source.df$cond.pi0.x, 
#                                                       p = p.vec, 
#                                                       pi0 = default1.source.df$est.pi0.x))
# 
# # default setting alpha = 1
# default2.source.df = Local.pi0.estimate.0127(x.vec, p.vec, neigh.m, cond.alpha = 1, B = 20)
# default2.q.values = Proposed.Q_0929(data = data.frame(orig.id = 1:m, 
#                                                       f_x = default2.source.df$cond.pi0.x, 
#                                                       p = p.vec, 
#                                                       pi0 = default2.source.df$est.pi0.x))
# 
# 
# sum(cv.q.values <= 0.2)
# sum(default1.q.values <= 0.2)
# sum(default2.q.values <= 0.2)

