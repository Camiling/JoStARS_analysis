rm(list=ls())
load('data/PanCancer_data.RData')
library(stabJGL)
library(foreach)

run_stabJGL = FALSE
run_FGL = TRUE
nCores = 20


# Use jointGHS on data ---------------------------------------------------------

if(run_stabJGL){
  set.seed(123)
  y = list(brca_dat, ovac_dat, ucec_dat)
  #res.joint = stabJGL::stabJGL(y, scale=T, var.thresh = 0.1, ebic.gamma = 1,lambda2.min=0, lambda2.max = 0.01, nCores=nCores, parallelize = T) 
  #res.joint = stabJGL::stabJGL(y, scale=T, var.thresh = 0.1, ebic.gamma = 1,lambda2.max = 0.05, nCores=nCores, parallelize = T) # best so far
  res.joint = stabJGL::stabJGL(y, scale=T, var.thresh = 0.1, ebic.gamma = 1,lambda2.max = 0.1,nlambda2=60, nCores=nCores, parallelize = T) 
  save(res.joint,file='data/res_stabJGL_PanCancer.RData')
}

if(run_FGL){
  source('useful_functions.R')
  set.seed(123)
  y = list(brca_dat, ovac_dat, ucec_dat)
  y = lapply(y,scale)
  res.fgl = JGL_select_AIC(Y=y,penalty='fused',nlambda1=20,lambda1.min=0.01,
                           lambda1.max=1,nlambda2=60,lambda2.min=0,lambda2.max=0.01,lambda2.init=0.01,penalize.diagonal=F)
  save(res.fgl,file='data/res_FGL_PanCancer.RData')
}

