#########################################################################################################
#### This R code implements the GSPPCA algorithm for high-dimensional unsupervised feature selection ####
#########################################################################################################

# References: 
# [1] P.-A. Mattei, C. Bouveyron and P. Latouche, Globally Sparse Probabilistic PCA, Proc. AISTATS 2016, pp. 976-984
# [2] C. Bouveyron, P. Latouche and P.-A. Mattei, Bayesian Variable Selection for Globally Sparse Probabilistic PCA, HAL preprint 01310409

# IMPORTANT REMARK: we use the model described in [2] rather than [1]. These models simply differ by the parametrisation of alpha.

# Contact: pierrealexandre.mattei@parisdescartes.fr 
# http://pamattei.github.io

library(MASS)
library(rARPACK)
library(Bessel)

tr=function(A) sum(diag(A))
trprod=function(A,B) sum(A*t(B)) # efficiently computes tr(A %*% B)
logdet=function(A) determinant(A,logarithm = TRUE)$modulus[1]

# This function implements the VEM algorithm described in section 2 of [2]
# 
# INPUT:
# data matrix x in R^{n*p} (each row is one observation)
# dimension of the latent space d in N*
# number of iterations nit in N*
# convergence criterion epsi>0 (the VEM stops when the relative difference of free energy becomes lower than epsi)
# initialization value alpha>0 (smaller alpha might lead to sparser solutions)
# initialization value s>0
# if verbose=TRUE, displays the iterations
# initialization values for u, Sigma, mu, S and M can also be specified (but it is usually unnecessary)
# 
# OUTPUT:
# relaxed parameter u in  [0,1]^p at convergence
# free energy monitored at each iteration L in R^nit
# value of the relative difference when the algorithm stops reldiff in R
# values of alpha,s,Sigma, mu, S and M at convergence


VEM<-function(x,d,nit=100,epsi=0.0001,alpha=1,s=1,verbose=FALSE,u=NULL,Sigma=NULL,mu=NULL,S=NULL,M=NULL){
  
  ## Data preprocessing
  x = scale(as.matrix(x),scale=FALSE)
  trcrossprodx = trprod(x,t(x))
  n = nrow(x)
  p = ncol(x)
  
  ## Hyperparameter initialisation
  if (is.null(u)) u = rep(1,p)
  if (is.null(alpha)) alpha = 1
  
  ## A posteriori parameters initialisation
  if (is.null(Sigma)) Sigma = diag(d)
  if (is.null(mu)|is.null(M)) dvs = svds(x,k=d)
  if (is.null(mu)) mu = dvs$u; #mu=matrix(rnorm(n*d),nrow=n)
  if (is.null(S)) S = array(rep(alpha^(-2)*diag(d),p),dim=c(d,d,p))
  if (is.null(M)) M = dvs$v; #M=matrix(rnorm(p*d),nrow=p)
  Su = array(do.call('c', lapply(1:p, function(i) u[i]^2*S[,,i])),dim=c(d,d,p))
  
  ## EM Algorithm
  L = numeric(nit) # the negative variational free energy
  j = 0
  reldiff = epsi+10
  
  while ((j<nit)&(reldiff>epsi)) {
    
    j=j+1
    
    ## E-Step
    Sigma = chol2inv(chol(diag(d)+s^(-2)*crossprod(u*M)+s^(-2)*apply(Su,c(1,2),sum)))
    
    mu= s^(-2)* x %*% (u * M) %*%t(Sigma)
    
    crossprodmu = crossprod(mu)
    
    S = array(do.call('c',lapply(1:p, FUN=function(k) chol2inv(chol( alpha^(2)*diag(d)+u[k]^2*s^(-2)*crossprodmu+u[k]^2*s^(-2)*n*Sigma ) ))),dim=c(d,d,p))
    
    M = t(sapply(1:p, FUN=function(k) s^(-2)*u[k]*S[,,k]%*%(t(mu)%*%x[,k])))
    if (d==1) M = t(t(sapply(1:p, FUN=function(k) s^(-2)*u[k]*S[,,k]%*%(t(mu)%*%x[,k]))))
    conc = do.call('c', lapply(1:p, function(i) u[i]^2*S[,,i]))
    Su = array(conc,dim=c(d,d,p))
    
    ## M-Step
    s2 = (n*p)^(-1)*(trcrossprodx+trprod(n*Sigma+crossprodmu,crossprod(u*M)+apply(Su,c(1,2),sum))-2*trprod(x%*%(u * M), t(mu)))
    if (d==1) s2 = (n*p)^(-1)*(trcrossprodx+trprod(n*Sigma+crossprodmu,crossprod(u*M)+apply(Su,c(1,2),sum))-2*trprod(x%*%(u * M), t(mu)))
    s = sqrt(s2)
    alpha = ((trprod(t(M),M) +tr(apply(S,c(1,2),sum)))/(d*p))^(-1/2)
    u = sapply(1:p, function(k){
      trprod(S[,,k]+tcrossprod(M[k,]),n*Sigma+crossprodmu)^(-1)*crossprod(M[k,],t(mu)%*%x[,k])
    })
    id0 = which(u<0)
    if(length(id0)>0) u[id0] = .Machine$double.eps
    id1 = which(u>1)
    if(length(id1)>0) u[id1] = 1
    
    ## Compute the free energy and check convergence
    L[j]=-n*p*log(s)-1/(2*s^2)*trcrossprodx + (d*p)*log(alpha) + s^{-2}*trprod(x, (u * M) %*% t(mu))-1/2*s^{-2}*trprod(n*Sigma+crossprodmu,crossprod(u*M)+apply(Su,c(1,2),sum)) -n/2*tr(Sigma)-1/2*tr(crossprodmu)-1/2*alpha^(2)*(trprod(M,t(M))+tr(apply(S,c(1,2),sum))) + n/2*logdet(Sigma) +1/2*sum(apply(S,c(3),logdet))
    
    if (verbose==TRUE) {cat(j); cat(' ')}
    if (j>10) reldiff=abs((L[j]-L[j-1])/L[j])
  }
  return(list(u=u,L=L[1:j],reldiff=reldiff,alpha=alpha,s=s,Sigma=Sigma,mu=mu,S=S,M=M))
}


# This function implements the model selection procedure described in section 2 of [2]
# 
# INPUT:
# data matrix x in R^{n*p} (each row is one observation)
# vector of variables weights u in [0,1]^p (in the case of the GSPPCA algorithm, u is the first output of the VEM, but other orderings can be tested)
# dimension of the latent space d in N*
# pmin and pmax are the minimum and maximum cardinalities of tested models
# method is the method used for noise variance estimation, can either be "JL" (Johnstone & Lu), "BC" (bias-corrected) or "ML" (maximum-likelihood)
#
#
# OUTPUT:
# qhat in [(d+1),p], number of relevant variables selected by the algorithm
# l, vector of variable length, value of the noiseless marginal likelihood for each model tested



Occam <- function(x,u,d,pmin=d+1,pmax=ncol(x)-1,method="BC",verbose=FALSE){
  
  ## Reorganize the variables according to the ordering of the coefficients of u
  xord <- scale(x[,order(u,decreasing=T)],scale=FALSE)
  n = nrow(x)
  p = ncol(x)
  
  ## Compute the estimator for the noise variance
  
  if ((method=="ML")|(method=="BC")) s = sdPPCA(xord,d,method=method)
  if (method=="JL") s=sqrt(median(apply(xord,MARGIN = 2,FUN=function (x0) sum(x0^2)/n)))
  
  ## Discard the variables whose u coefficient is exactly zero (ARD effect)
  
  pmax = min(pmax,p-sum(u<.Machine$double.eps))
  
  dims=seq(pmin,pmax,1) # cardinalities of the models tested
  
  ## Define a function to compute the marginal likelihood of a model of cardinality k
  
  computelogmarginal<- function (k) {
    if (k==1){ R=sqrt(xord[,1]^2) }else{ R=sqrt(rowSums(xord[,1:k]^2)) }
    
    # Define the function to optimize and its derivative
    neglogl<- function(x0){
      if((k-d)>10){logevi = sum(besselK.nuAsym(R*x0,(k-d)/2,log=TRUE,k.max=4,expon.scaled = TRUE)) -sum(R)*x0  +(d-k)/2*sum(log(R)) + n*(-(k+d)/2+1)*log(2) - n*lgamma(d/2) - n*k/2*log(pi) + n*(k+d)/2*log(x0) -1/2*sum(xord[,(k+1):p]^2)/s^2  - n*(p-k)*log(sqrt(2*pi)*s)}
      else {logevi = sum(log(BesselK(R*x0,(k-d)/2,expon.scaled = TRUE)))-sum(R)*x0  +(d-k)/2*sum(log(R)) + n*(-(k+d)/2+1)*log(2) - n*lgamma(d/2) - n*k/2*log(pi) + n*(k+d)/2*log(x0) -1/2*sum(xord[,(k+1):p]^2)/s^2  - n*(p-k)*log(sqrt(2*pi)*s)}
      return(-logevi)
    }
    
    Dneglogl<- function(x0){
      if((k-d)>10){Dlogevi = sum(-R*exp(besselK.nuAsym(R*x0,(k-d)/2-1,log=TRUE,k.max=4,expon.scaled = TRUE)-besselK.nuAsym(R*x0,(k-d)/2,log=TRUE,k.max=4,expon.scaled = TRUE))) +n*d/x0}
      else {Dlogevi = sum(-R*BesselK(R*x0,(k-d)/2-1,expon.scaled = TRUE)/(BesselK(R*x0,(k-d)/2,expon.scaled = TRUE))) + n*d/x0}
      return(-Dlogevi)
    }
    return(-optim(0.1,fn=neglogl,gr=Dneglogl,method ="L-BFGS-B",lower = 10^(-10), upper = Inf)$value) }
  
  ## Compute the noiseless marginal likelihood  
  
  l = sapply(dims,computelogmarginal)
  if (verbose==TRUE) l = sapply(dims,computelogmarginal)
  
  
  return(list(qhat=dims[which.max(l)],l=l))
}

# This function implements the whole GSPPCA algorithm
#
# INPUT:
# data matrix x
# dimension of the latent space d
# some other parameters from the VEM and Occam functions (in particular, we can choose several starting points for the alpha parameter, each one is used to run 3 iterations and the best one is chosen)
# numvar="Occam" indicates that you use the noiseless marginal for model selection
# but you can also use an integer argument if you already know the number of relevant variables
#
# OUTPUT:
# sparsity pattern vhat in {0,1}^p
# sparsity level qhat (number of relevant variables)
# if numvar="Occam", the values of the noiseless marginal likelihood for the tested models
# loadings and scores are globally sparse PCA loadings and scores obtained after performing PCA on the relevant variables
# center is the centering used before performing PCA


gsppca <- function (x,d,nit=200,alpha=c(1,10),s=1,epsi=1e-5,numvar="Occam",pmin=(d+1),pmax=ncol(x)-1,method="BC",verbose=FALSE){
  
  n <- nrow(x); p <- ncol(x)
  
  # Try several initializations values for alpha and take the best one
  if (length(alpha)>1){
    L <- numeric(length(alpha))
    numzeros <- numeric(length(alpha))
    outs <- lapply(alpha, function(a) VEM(x=x,d=d,nit=3,alpha=a,s=s,epsi=epsi,verbose=FALSE))
    for (i in  1:length(alpha)) {
      L[i] <- outs[[i]]$L[length(outs[[i]]$L)]
      numzeros[i] <- sum(outs[[i]]$u<.Machine$double.eps)
    }
    if (sum(numzeros<=(p-d))==0) stop("All initial values for alpha led to u=0 ! Choose different starting points.")
    bestout <- outs[[which.max(L[numzeros<=(p-d)])]]
    rm(outs)
  }
  
  
  
  # Run the VEM
  if (verbose==TRUE) cat('Iterations ')
  if (length(alpha)>1) out <- VEM(x=x,d=d,nit=nit,alpha=bestout$alpha,s=bestout$s,epsi=epsi,verbose=verbose,u=bestout$u,Sigma=bestout$Sigma,mu=bestout$mu,S=bestout$S,M=bestout$M)
  if (length(alpha)==1) out <- VEM(x=x,d=d,nit=nit,alpha=alpha,verbose=verbose)
  
  u <- out$u
  
  if (sum(u<.Machine$double.eps)>=(p-d)) stop("VEM converged to u=0 !")
  
  # Model selection
  if (numvar=="Occam") { if (verbose==TRUE) {cat('\n'); cat('Model Selection...\n')}
    res <- Occam(x=x,u=u,d=d,pmin=pmin,pmax=pmax,method=method,verbose=verbose)
    if (verbose==TRUE) cat('done.')
    qhat<-res$qhat
    l<-res$l}
  else qhat<-numvar
  vhat = rep(0,p)
  vhat[order(u,decreasing=TRUE)[1:qhat]] = 1
  
  # PCA on the selected variables
  outpca = prcomp(x[,vhat==1],retx = )
  loadings = matrix(0,ncol=d,nrow=p)
  loadings[vhat==1,]<-outpca$rotation[1:d]
  center = colMeans(x)
  
  if (numvar=="Occam") {
    return(list(l=l,vhat=vhat,qhat=qhat,u=u,loadings=loadings,scores=outpca$x[,1:10],center=center))}
  else return(list(vhat=vhat,qhat=qhat,loadings=loadings,scores=outpca$x[,1:10],center=center))
}



# This function computes an estimator of the standard deviation of the noise in a probabilistic PCA model
#
# REFERENCES :
#
# [1] D. Passemier, Z. Li, J. Yao, On estimation of the noise variance in high-dimensional probabilistic principal component analysis, JRSSB 2016
# [2] M. Tipping and C. Bishop, Probabilistic principal component analysis, JRSSB 1999
# 
# INPUT:
# data matrix x in R^{n*p} (each row is one observation)
# dimension of the latent space d in N*
# method can be either "BC" which corresponds to the bias-corrected estimator of [1] or "ML" which corresponds to the classical ML estimator of [2]
#
# OUTPUT:
# noise standard error estimator s
#
# Contact: pierrealexandre.mattei@parisdescartes.fr 
# http://pamattei.github.io/

sdPPCA <- function(X,d,method="BC"){
  x <- scale(X,scale=FALSE)
  
  n <- as.numeric(nrow(x))
  p <- as.numeric(ncol(x))
  
  ratio <- as.numeric(p/n)
  
  E <- c(rep(0,p-min(n,p)),sort(svd(x,nu=0,nv=0)$d^2)/n)
  sML2 <- sum(E[1:(p-d)])/(p-d)
  
  alpha=numeric(d)
  for (j in 1:d){
    if((sML2*(ratio+1)-E[p-j+1])^2-4*ratio*sML2^2<0)  stop("BC method for noise estimation failed: no real root in the polynomial! Try another method like JL or ML.")
    alpha[j] = 0.5*(sqrt((sML2*(ratio+1)-E[p-j+1])^2-4*ratio*sML2^2)-sML2*(ratio+1)+E[p-j+1])}
  b=sqrt(ratio/2)*(d+sML2*sum(1/alpha))
  
  sPY2=sML2+(b*sML2*sqrt(2*ratio))/(p-d)
  
  if (method=="BC") s=sqrt(sPY2)
  if (method=="ML") s=sqrt(sML2)
  
  return(s=s)
  
  
}

