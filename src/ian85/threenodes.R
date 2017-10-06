# Look at marginal pairwise distributions over fully connected graphs with given graph weights. 

# These weights are probabilities for, in order, 1 = 2, 1 = 3, and 2 = 3. Feel free to tinker.
weight.vec = c(.9, .1, .9)

# I was investigating the eigendecomposition of the 'transition matrix', if you like, but didn't really find anything, so the matrix isn't useful atm.
weight.matrix = diag(3)
weight.matrix[upper.tri(weight.matrix)] = weight.matrix[lower.tri(weight.matrix)] = weight.vec
eigen(weight.matrix)

# The full pmf is a distribution over the three nodes having the labels 1, 2, or 3, so 27 total combinations (some isomorphic, though). Interestingly, the number of non-isomorphic labelings is the third Bell number, which is a fun curio but probably not useful: https://en.wikipedia.org/wiki/Bell_number

# This makes a pmf table and computes the probabilities for each configuration
pmf = as.matrix(expand.grid(1:3, 1:3, 1:3))
probs = numeric(nrow(pmf))
for(i in 1:nrow(pmf)){
  bool = c(pmf[i,1] == pmf[i,2], pmf[i,1] == pmf[i,3], pmf[i,2] == pmf[i,3])
  probs[i] = prod(abs(weight.vec - !bool))
}

# This is some mcmc code to sample from the full joint distribution for the network configurations. Each row of mcmc is a sample. Empirical summaries should match the theoretical values in probs, which they seem to do

# number of mcmc runs
nsamp = 100000

# This randomly pulls rows from the pmf
mcmc = pmf[sample(1:nrow(pmf), nsamp, prob = probs, rep = T), ]
names(mcmc) = NULL

# This counts the number of times 1 = 2, 1 = 3, and 2 = 3.
same.props = c(mean(mcmc[,1] == mcmc[,2]), mean(mcmc[,1] == mcmc[,3]), mean(mcmc[,2] == mcmc[,3]))
same.props

# Counts the number of times c(1, 1, 1) appears in the sample. Should be close to probs[1] AFTER NORMALIZATION
mean(apply(mcmc, 1, function(x) identical(as.numeric(x), c(1, 1, 1))))
all.equal.prob = mean(apply(mcmc, 1, function(x) identical(as.numeric(x), c(1, 1, 1))))
all.equal.prob
probs[1]/sum(probs)

# TO DO:
# write a function to compute the full empirical pmf and compare to the actual pmf.
