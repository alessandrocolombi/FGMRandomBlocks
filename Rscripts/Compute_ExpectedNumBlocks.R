theta = 2
sigma = 0.0001
p = 40

A = gamma(theta)/gamma(theta+sigma) * gamma(theta+p+sigma)/gamma(theta+p) - 1

logK = log(theta) - log(sigma) + log(A)
K = exp(logK)
K
