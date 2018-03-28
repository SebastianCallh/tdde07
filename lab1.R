# 1. Bernoulli ... again.

# a) Show that posterior converges to 14/20 approx. 0.7
alpha0 <- 2
beta0  <- 2 
nDraws <- 10000
n      <- 20
s      <- 14
f      <- n - s
post   <- rbeta(nDraws, alpha0 + s, beta0 + f)
hist(post, breaks = 100)

# b) Simulation to compute posterior probability P(theta < 0.4 | y)
thetas <- rbeta(nDraws, alpha0 +s , beta0 + f)
length(thetas[thetas<0.4]) / length(thetas)  # Simulated value
pbeta(0.4, alpha0+s, beta0+f)                # Exact value
hist(thetas, breaks = 100)

# c) Log odds posterior distribution
fi <- sapply(thetas, function(theta) { log(theta / (1 - theta)) })
hist(fi, breaks = 100)

# 2. Log-normal distribution and the Gini coefficient.
incomes <- c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
mu <- 3.5
n  <- length(incomes)
tau2 <- sum((log(incomes) - mu)^2) / n 

# a) Simulate posterior draws
X      <- rchisq(nDraws, n)
sigma2 <- n*tau2/X
hist(sigma2, breaks = 100)

# Define chi squared in terms of Gamma distribution
dinvgamma <- function(x, a, b) {
  (b^a)/gamma(a) * x^(-a-1) * exp(-b/x)
}

dinvchisq2 <- function(x, n, tau2) {
  a <- n / 2
  b <- n*tau2 / 2
  dinvgamma(x, a, b)
}

plot(dinvchisq2(seq(0, 3, 0.01), n, tau2), type = 'l') # Distribution

# b) Gini coefficient posterior
G      <- 2*pnorm(sqrt(sigma2/2)) - 1
G.hist <- hist(G, breaks = 1000, xlim = 0:1)

# c) G Credible interval
G.slice <- sort(G)[(0.025*nDraws):(0.975*nDraws)]
cred.lower.bound <- min(G.slice)
cred.upper.bound <- max(G.slice)

# c) G Highest Posterior Density (HPD)
G.dens  <- density(G)
plot(G.dens, type = 'p')

i <- order(G.dens$y, decreasing = TRUE)
csum <- cumsum(G.dens$y[i])
m    <- length(csum[csum < sum(G.dens$y)*0.95])
H    <- G.dens$x[i][1:m]
hpd.lower.bound <- min(H)
hpd.upper.bound <- max(H)

# 3
# Each row in matrix p ~ 1 value of k, wrong, should be scalar
Y  <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
K  <- dexp(seq(0, 5, 0.01))
mu <- 2.39


Von.Mises <- function (y, mu, k){
  exp(k*cos(y - mu)) / (2*pi*besselI(k, 0))
}

likelihood <- function (k) { prod(Von.Mises(Y, mu, k)) }
posterior  <- function (k) { likelihood(k) * k } #?

plot(K)
plot(sapply(K, likelihood))
plot(sapply(K, posterior))


plot(Von.Mises(Y, mu, 0), ylim = c(0,4))
points(Y)


# Sebas kvällsmeck
Y          <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mu         <- 2.39
lambda     <- 1
kappas     <- seq(0, 10, 0.05)
likelihood <- sapply(kappas, function (k) { prod(Von.Mises(Y, mu, k)) })
prior      <- dexp(kappas, lambda)
posterior  <- likelihood * prior
plot(posterior, col = 'orange', type = 'l')
kappas[which.max(posterior)]
#plot(kappas, type ='l', col = 'green')
#lines(sapply(kappas, likelihood), col = 'blue')
