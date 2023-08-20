estimate.m0=function(p, B = 20)
{
  #
  #This function estimates the number of true null hypotheses given a vector of p-values
  #using the method of Nettleton et al. (2006) JABES 11, 337-356.
  #The estimate obtained is identical to the estimate obtained by the iterative
  #procedure described by Mosig et al. Genetics 157:1683-1698.
  #The number of p-values falling into B equally sized bins are counted.
  #The count of each bin is compared to the average of all the bin counts associated
  #with the current bins and all bins to its right.  Working from left to right, 
  #the first bin that has a count less than or equal to the average is identified.
  #That average is multiplied by the total number of bins to obtain an estimate of m0, 
  #the number of tests for which the null hypothesis is true.
  #
  m <- length(p)
  m0 <- m
  bin <- c(-0.1, (1:B)/B)
  bin.counts=rep(0,B)
  for(i in 1:B){
    bin.counts[i]=sum((p>bin[i])&(p<=bin[i+1]))
  }
  tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
  temp <- bin.counts - tail.means
  index <- min((1:B)[temp <= 0])
  m0 <- B * tail.means[index]
  return(m0)
}

jabes.q = function(p,B=20)
{
  #
  #This function computes q-values using the approach of Nettleton et al.
  #(2006) JABES 11, 337-356.
  #
  #Author: Dan Nettleton
  #
  
  m = length(p)
  m0=estimate.m0(p,B)
  k = 1:m
  ord = order(p)
  p[ord] = (p[ord] * m0)/(1:m)
  qval = p
  qval[ord]=rev(cummin(rev(qval[ord])))
  
  return(qval)
}


fdr_0927 = function(p,B=20)
{
  #
  #This function computes q-values using the approach of Nettleton et al.
  #(2006) JABES 11, 337-356.
  #
  #Author: Dan Nettleton
  #
  
  m = length(p)
  m0=estimate.m0(p,B)
  k = 1:m
  ord = order(p)
  p[ord] = (p[ord] * m0)/(1:m)
  qval = p
  # qval[ord]=rev(cummin(rev(qval[ord])))
  
  return(qval)
}

######################################################################################

orr.estimate.m0 = function(p1, p2, B = 20){
  
  # Settings
  # p = data.frame(list_p %>% do.call("cbind",.))
  # colnames(p) = paste0("p",1:ncol(p))
  
  p = cbind(p1, p2)
  m = nrow(p)
  k = ncol(p)
  bin = seq(0, 1, 1/B)
  # bin.counts and tail.means
  bin.info = apply(t(sapply(bin, function(b){
    counts = apply(p,2,function(.){sum(.<=b)})
    return(c(counts))
  })), 2, diff)
  tail.info = apply(bin.info, 2, function(.){rev(cumsum(rev(.))/(1:B))}) 
  
  indices = sapply(1:k, function(i){
    min(which(bin.info[,i]<=tail.info[,i]))
  })
  
  cutoffs = bin[-(B+1)][indices] # V treating zero...?
  # print(cutoffs)
  m0 =  sum(eval(parse(text = paste(paste0("(",colnames(p),">",cutoffs,")"),collapse="&"))))/prod(1-cutoffs)
  
  return(m0)

}

mod.jabes.q = function(p, pi0, B = 20) {

  m = length(p)
  m0 = m*pi0
  k = 1:m
  ord = order(p)
  p[ord] = (p[ord] * m0)/(1:m)
  qval = p
  qval[ord] = rev(cummin(rev(qval[ord])))
  
  return(qval)
}

######################################################################################

