# 1. Bernoulli ... again.
library(geoR)
# a) Show that posterior converges to 14/20 approx. 0.7
alpha0 <- 2
beta0  <- 2 
nDraws <- 10000
n      <- 20
s      <- 14
f      <- n - s
post   <- rbeta(nDraws, alpha0 + s, beta0 + f)
hist(post)

# b) Simulation to compute posterior probability P(theta < 0.4 | y)
thetas <- rbeta(nDraws, alpha0 +s , beta0 + f)
length(thetas[thetas<0.4]) / length(thetas)  # Simulated value
pbeta(0.4, alpha0+s, beta0+f)                # Exact value
hist(thetas)

# c) Log odds posterior distribution
fi <- sapply(thetas, function(theta) { log(theta / (1 - theta)) })
hist(fi)  

# Ask question: Did we the fuck do this good? Or right? or somethin?



# 2. Log-normal distribution and the Gini coefficient.
incomes <- c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
mu <- 3.5
n  <- length(incomes)
tau2 <- sum((log(incomes) - mu)^2) / n 

# a) simulate draws
# Do the histograms need normalization before any conclusions?
draws <- rinvchisq(nDraws, n, tau2)
x     <- seq(0, 3, 0.01)
hist(draws, breaks = 1000)  # Simulated 
plot(dinvchisq(x, n, tau2)) # Distribution

# b) Gini coefficient posterior
# If one bin contains a value such that the cumsum exceeds i.e. 25 it will still be dropped
# What does 1st Qu.3rd Qu.do?
G      <- 2*pnorm(sqrt(draws/2)) - 1
G.hist <- hist(G, breaks = 1000, xlim = 0:1)
csum   <- cumsum(G.hist$counts)
start  <- length(csum[csum < 250])
end    <- length(csum[csum < 9750])
shit <- G.hist$counts[start:end]
plot(G.hist$counts[start:end], type = "h")

# c) G credibility interval
H <- G / sum(G)
H.hist <- hist(H, breaks = 1000)

csum <- cumsum(sort(G / sum(G)))
hist(sort(G / sum(G)))
start <- length(csum[csum < 0.025])
end   <- length(csum[csum < 0.975])

G_prime <- G
G_prime <- G_prime[-(1:start)]
hist(G_prime, breaks = 1000, xlim = 0:1)
#G_prime <- G_prime[-(end:length(G))]

density(G)

# 3
# Each row in matrix p ~ 1 k value
y  <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
k  <- dexp(seq(0, 5, 0.01))
mu <- 2.39
Von.Mises <- function (y, mu, k){
  exp(k*cos(y - mu)) / (2*pi*besselI(k, 0))
}
p <- sapply(k, function (x) { Von.Mises(y, mu, x) })
shist <- hist(p, breaks = 100)
density(p)
