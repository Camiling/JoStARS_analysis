library(foreach)
library(stats)
library(ggplot2)
library(GGally)
library(sna)
library(JGL)
library(gridExtra)
library(igraph)
library(network)
library(superheat)
library(ggfortify)
library(doParallel)
library(foreach)
library(parallel)
library(dplyr)

# Find Matthews correlation coefficient for estimated graph
MCC = function(g,g.hat){
  p = nrow(g[,])
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat) = rep(0,p) 
  tp = sum(g.hat ==1 & g ==1)/10 # True positives. Divide by 10 to avoid integer overflow. 
  fp = sum(g.hat ==1 & g ==0)/10 # False positives
  tn = (sum(g.hat == 0 & g == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat == 0 & g == 1)/10 # False negatives
  return((tp*tn - fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
}


get_and_print_edges <- function(a.mat, col.names, theta.mat=NULL) {
  # Function for printing all the edges in a graph with a layout that can be inserted into a latex table.
  # Also returns a data frame containing the edges
  # a.mat:          the adjacency matrix
  # col.names:      the names of the nodes in the graph
  # theta.mat:      the precision matrix. If included, the size the partial correlations are included in the table as well
  a.mat[which(diag(rep(1, ncol(a.mat))) == 1, arr.ind = T)] <- 0 # Make diagonal zero
  pairs <- which(a.mat[, ] == 1, arr.ind = T)
  df <- data.frame(t(apply(pairs, 1, sort))) # Sort so that the node in the pair whose name is first in the alphabet is first.
  df <- unique(df)
  names <- cbind(col.names[df[, 1]], col.names[df[, 2]])
  if (!is.null(theta.mat)){
    effect = round(-cov2cor(theta.mat)[cbind(df[,1], df[,2])], 5)
    names = cbind(names, effect)
    return(names)
  }
  for (i in 1:nrow(names)) {
    cat(names[i, 1], " & ", names[i, 2], " \\\\ \n")
  }
  return(names)
}



JGL_select_AIC = function(Y,penalty='fused',nlambda1,lambda1.min,lambda1.max,nlambda2,lambda2.min,lambda2.max, 
                          lambda2.init,penalize.diagonal){
  # JGL with parameters selecetd by the adapted AIC crierion (Danaher et al.)
  K=length(Y)
  p=ncol(Y[[1]])
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,cov)
  lambda1.vals = seq(lambda1.min,lambda1.max,length.out=nlambda1)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out=nlambda2)
  mods.lam1=list()
  aic.lam1=rep(0,length(lambda1.vals)) 
  for (i in 1:length(lambda1.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1.vals[i],lambda2 = lambda2.init,return.whole.theta = T,penalize.diagonal=penalize.diagonal)$theta
    mods.lam1[[i]] = mod.temp
    aic.lam1[i] = AIC_adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind.lam1 = which.min(aic.lam1)
  lambda1=lambda1.vals[opt.ind.lam1]
  
  mods.lam2=list()
  aic.lam2=rep(0,length(lambda2.vals))
  for (i in 1:length(lambda2.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1,lambda2 = lambda2.vals[i],penalize.diagonal = penalize.diagonal,
                   return.whole.theta = T)$theta
    mods.lam2[[i]] = mod.temp
    aic.lam2[i] = AIC_adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind = which.min(aic.lam2) 
  res=list(opt.fit=mods.lam2[[opt.ind]],opt.lambda1 = lambda1,opt.lambda2 = lambda2.vals[opt.ind],
           opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
  return(res)
}





