set.seed(123)

# First, we simulate some data according to a simple globally sparse PPCA model

n = 100
p = 500
d = 10 # dimension of the latent space
q = 50 # true number of relevant variables

vtrue = c(rep(0,p-q),rep(1,q)) # true sparsity pattern
sigma = 1
W = diag(vtrue) %*% matrix(rnorm(p*d),nrow=p)
X = matrix(rnorm(n*d),nrow=n) %*% t(W) + sigma*matrix(rnorm(p*n),nrow=n)


# Now we can run the GSPPCA algorithm
# On this 100x500 example, it takes less than 10 seconds on a MacBook Pro i5 processor

out = gsppca(X,d=10,epsi=1e-4,method="JL",verbose=TRUE)

# and check if the true sparsity pattern is recovered

print(sum(abs(out$vhat - vtrue)) == 0)

# One can plot the values of the marginal likelihood for growing cardinalities

plot((d+1):(d+length(out$l)),out$l,xlab='Cardinality',ylab='Noiseless evidence',type='l')

# As a simple illustration, we can compare it with SPCA (Zou, Hastie & Tibshirani, JCGS 2006)

library(elasticnet)
vhatspca = spca(X,K = 1,sparse = "varnum",para = q)$loadings!=0

# We can see that even if we choose the true number of relevant variables, SPCA makes many mistakes

print(sum(abs(vtrue-vhatspca)))
