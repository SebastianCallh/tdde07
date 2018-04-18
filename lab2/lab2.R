library(mvtnorm)
data <- read.csv("TempLinkoping.txt", sep = "\t")
X <- cbind(rep(1, length(data$time)), data$time, data$time^2)
Y <- data$temp

# 1. Linear and polynomial regression

rinvchisq <- function(vx, sigmax, draws) {
  vx*sigmax/rchisq(draws, vx)
}

#a) Determining the prior distribution
delta <- 0.01
times <- seq(delta, 1, delta)
nDraws  <- 5000 
posteriorDraws <- 5000
mu0    <- as.matrix(c(-10, 150, -150))   # From tuning the plotted polynomial
omega0 <- 1*diag(3)           # same variance for betas, we dont know much
v0      <- 20
sigma0  <- 2
sigma2s <- rinvchisq(v0, sigma0, nDraws)#v0*sigma0/rchisq(nDraws, v0)
n       <- dim(X)[1]

glmthing <- lm(temp ~ time + I(time^2), data = data)
summary(glmthing)

# Plot of prior variance
plot.prior.varaince <- function() {
  plot(density(sigma2s),
       main = "Sigma^2 prior",
       xlab = expression(sigma^2),
       ylab = "Density",
       type = 'l')  
}

# Thetas are matrix of theta - all parameters with priors
thetas <- cbind(t(sapply(sigma2s, function(sigma2) { 
  rmvnorm(1, mu0, sigma2*solve(omega0))
})), sigma2s)

# y is the prior distribution of the model
prior.temp <- apply(thetas[seq(1, 10),], 1, function(theta) {
  sapply(times, function(time) {
    theta[1] + theta[2]*time + theta[3]*time^2 #+ rnorm(1, 0, theta[4])
  })
})


#b) Check if prior distribution is sensible
m <- dim(prior.temp)[1]
x.axis <- (1:m) / m

plot.prior.regression.curves <- function () {
  plot(x.axis, prior.temp[,1], type = 'l', ylim = c(-15, 40))
  for (i in 2:10) {
    lines(x.axis, prior.temp[,i])
  }
}

# c) Simulating from posterior distribution

# Posterior mapping
beta.hat <- solve(t(X)%*%X)%*%t(X)%*%Y # OLS 
A   <- t(X) %*% X
B   <- as.numeric(t(Y) %*% Y)
mun      <- solve(A + omega0) %*% (A%*%beta.hat + omega0%*%mu0)
omegan   <- A + omega0
vn       <- v0 + n
sigma2n  <- as.numeric((v0%*%sigma0 + (B + t(mu0)%*%omega0%*%mu0 - t(mun)%*%omegan%*%mun))/vn)

post.sigma2s <- rinvchisq(vn, sigma2n, posteriorDraws)
post.thetas  <- cbind(t(sapply(post.sigma2s, function(post.sigma2) {
  rmvnorm(1, mun, post.sigma2 * solve(omegan)) 
})), post.sigma2s)

post.betas <- post.thetas[,1:3]
theta.mean <- apply(post.thetas, 2, mean)
post.pred.mean <- X  %*% as.matrix(theta.mean[1:3])

nSamples <- 5000
bounds <- as.matrix(sapply(times, function(t) {
  ys <- as.matrix(apply(post.betas[1:nSamples,], 1, function(beta) {
    c(1, t, t^2) %*% as.matrix(beta)
  }))
  
  lower.bound <- nSamples * 0.05
  upper.bound <- nSamples * 0.95
  bounded.y   <- sort(ys)[lower.bound:upper.bound]
  list(upper = head(bounded.y, 1), lower = tail(bounded.y, 1))
}))

plot.posterior.betas <- function() {
  plot(data$time, data$temp, ylim=c(-20,25))
  lines(data$time, post.pred.mean)
  lines(x=times,y=bounds[1,], col = 'Blue')
  lines(x=times,y=bounds[2,], col = 'Red')
}

# d)
x_tilde <- which.max(post.pred.mean)
# Hottest day of the year (# 197)


# Distribution over the hottest day
day <- as.matrix(apply(post.betas, 1, function(beta) {
  (-beta[2] / (2*beta[3]))*n
}))
mean(day)
plot.day.density <- function() {
  plot(density(day))
}

# e)
# mu0: set the last 4 values to 0 to delete the effect of the higher order terms
# omega: set the 4 last values to very large to prevent high variance for the higher order terms


# 2. Posterior approximation for classification with logistic regression

# a)
data2 <- read.table("WomenWork.dat", header=TRUE)
y    <- as.vector(data2[,1])
X    <- as.matrix(data2[,-1])

glmModel <- glm(Work ~0 + ., data=data2, family = binomial)
#summary(glmModel)

# b)
# want posterior
# to get that first we need prior and likelihood

beta.log.posterior <- function(betas, X, y, mu, sigma) {
  d    <- length(betas)
  pred <- X %*% betas
  
  log.likelihood <- sum(y*pnorm(pred, log.p = TRUE) + (1-y)*pnorm(pred, log.p = TRUE, lower.tail = FALSE))
  log.prior      <- dmvnorm(betas, matrix(mu, d, 1), sigma*diag(d), log = TRUE)
  
  log.likelihood + log.prior
}

tau2 <- 10^2
mu   <- 0
d <- dim(X)[2]
beta.initial <- rep(0, d)
result       <- optim(beta.initial, beta.log.posterior, gr = NULL, 
                      X, y, mu, tau2, 
                      hessian = TRUE,
                      method  = c("BFGS"),
                      control = list(fnscale=-1))

beta_hat <- result$par
J        <- result$hessian 

# numerical values should be in the report
beta.post <- rmvnorm(10000, beta_hat, -solve(J))

#NSmallChild
nsc <- beta.post[,7]
quantile(nsc, c(0.025, 0.975))
plot.density.nsc <- function() {
  plot(density(nsc))
}

# compute 95% credible interval for NSmallChild?
x <- c(1, 10, 8, 10, 1, 40, 1, 1)
y_hat <- x %*% t(beta.post)
p <- 1 / (1 + exp(-y_hat))
plot.does.work <- function() {
  plot(density(p))
  quantile(p, c(0.025, 0.975))
}

