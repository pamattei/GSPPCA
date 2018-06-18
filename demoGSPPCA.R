
source("GSPPCA.R")

set.seed(1234)

# First, we simulate some data according to a simple globally sparse PPCA model

n = 50
p = 500
d = 10 # dimension of the latent space
q = 20 # true number of relevant variables

vtrue = c(rep(0,p-q),rep(1,q)) # true sparsity pattern
sigma = 1
W = diag(vtrue) %*% (matrix(rnorm(p*d),nrow=p))
X = matrix(rnorm(n*d),nrow=n) %*% t(W) + sigma*matrix(rnorm(p*n),nrow=n)


# Now we can run the GSPPCA algorithm
# On this 100x500 example, it takes less than 10 seconds on a MacBook Pro i5 processor

out = gsppca(X,d=10,epsi=1e-4,method="BC",verbose=TRUE)

# and check if the true sparsity pattern is recovered

print(sum(abs(out$vhat - vtrue)) == 0)

# One can plot the values of the marginal likelihood for growing cardinalities

plot((d+1):(d+length(out$l)),out$l,xlab='Cardinality',ylab='Noiseless evidence',type='l')

# As a simple illustration, we can compare it with SPCA (Zou, Hastie & Tibshirani, JCGS 2006)

library(elasticnet)
vhatspca = spca(X,K = 1,sparse = "varnum",para = q)$loadings!=0

# We can see that even if we choose the true number of relevant variables, SPCA makes a few mistakes

print(sum(abs(vtrue-vhatspca)))


# Let's illutrate the advantages of the informative prior by generating very HD and very noisy data:

n = 50
p = 5000
d = 5 # dimension of the latent space
q = 20 # true number of relevant variables

vtrue = c(rep(0,p-q),rep(1,q)) # true sparsity pattern
sigma = 2
W = diag(vtrue) %*% (matrix(rnorm(p*d),nrow=p))
X = matrix(rnorm(n*d),nrow=n) %*% t(W) + sigma*matrix(rnorm(p*n),nrow=n)

# The uniform prior selects too many variables:

out = gsppca(X,d=d,epsi=1e-6,method="JL",verbose=TRUE)
plot((d+1):(d+length(out$l)),out$l,xlab='Cardinality',ylab='Noiseless evidence',type='l',main="Uniform prior over model space") # Evidence plot
plot((d+1):(d+200),out$l[1:200],xlab='Cardinality',ylab='Noiseless evidence',type='l') # with a zoom
print(sum(abs(out$vhat-vtrue)))


# In contrast, the sparsity-inducing prior works much better

outinform = gsppca(X,d=d,epsi=1e-6,nit = 1000,method="JL",prior="info",verbose=TRUE)
plot((d+1):(d+length(outinform$l)),outinform$l,xlab='Cardinality',ylab='Noiseless evidence',type='l',main="Sparsity-inducing prior over model space") # Evidence plot
plot((d+1):(d+200),outinform$l[1:200],xlab='Cardinality',ylab='Noiseless evidence',type='l') # with a zoom
print(sum(abs(outinform$vhat-vtrue)))


# Again, even with the right sparsity pattern, SPCA makes more mistakes than GSPPCA with the informative prior

vhatspca = spca(X,K = 1,sparse = "varnum",para = q)$loadings!=0
print(sum(abs(vtrue-vhatspca)))






