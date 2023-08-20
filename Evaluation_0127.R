library(pROC)

# # example code
# qvalue = runif(10000)
# true.degs = sort(which(runif(10000)<0.3))

pAUC_0127 = function(qvalue, true.degs){

m = length(qvalue)

# A vector of control-case status.
response = rep("non-DEG", m)
response[true.degs] = "true-DEG"

# Make the response as a factor with levels = c(control, case)
response = factor(response, levels = c("non-DEG","true-DEG"))

# If qvalue <= t, it is declared positive => direction = ">"
# FPR = c(0, 0.1) => specificity = c(1, 0.9) => partial.auc = c(1, 0.9)

T.roc = roc(response = response, predictor = qvalue, direction = ">")

# AUC value :
AUC = as.numeric(auc( roc = T.roc ))
# pAUC value :
pAUC = as.numeric(auc( roc = T.roc,
                       partial.auc = c(1, 0.9),
                       partial.auc.focus = "specificity",
                       partial.auc.correct = T))

# pAUC = as.numeric(auc( roc = roc(response = response, predictor = qvalue, direction = ">")))
return(c(AUC = AUC, pAUC = pAUC))

}

# Motivated from Evaluation_0928 # Almost the same !!

Evaluation_0127 = function(qvalue, true.degs, alpha){

  # IDs declared significant :
  positive.ids = which( qvalue <= alpha )

  # Evaluation based on three statistics :
  # True discovery, False discovery proportion, and pAUC
  
  TD = sum(positive.ids %in% true.degs)
  FDP = sum(!positive.ids %in% true.degs)/max(length(positive.ids),1)
  v.AUC = pAUC_0127(qvalue, true.degs)

  v.summary = c( TD = TD, FDP = FDP, v.AUC )

  return(v.summary)
  
}


