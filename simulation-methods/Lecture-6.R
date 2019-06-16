# Simulation Methods in Financial Engineering
# R ode for Lecture 6
# Brownian motion, geometric Brownian motion and their variants
# Jieli Shen


# 1. Standard Brownian motion
# Random walk construction
t <- seq(0, 500, 0.5)
n <- length(t) - 1
tmp <- matrix(t[-1] - t[-n], nrow = n, ncol = n, byrow = T)
tmp[row(tmp) < col(tmp)] <- 0
A <- sqrt(tmp)
z <- rnorm(n)
x <- A %*% z
x <- c(0, x)
plot(t, x, type = "l")
# Random walk alternative construction
t <- seq(0, 500, 0.5)
n <- length(t)
x <- rep(0, n)
for (i in 2:n) {
	x[i] <- x[i - 1] + sqrt(t[i] - t[i - 1]) * rnorm(1)
}
plot(t, x, type = "l")

# 2. Brownian motion generla form
mu <- 0.5
sigma <- 3
t <- seq(0, 500, 0.5)
n <- length(t)
x <- rep(0, n)
for (i in 2:n) {
	x[i] <- x[i - 1] + mu * (t[i] - t[i - 1]) + sigma * sqrt(t[i]- t[i - 1]) * rnorm(1)
}
plot(t, x, type = "l")

# 3. Brownian motion with time varying parameters
t <- seq(0, 500, 0.5)
n <- length(t)
mu <- sin(t / 50)
sigma <- sqrt(9 + cos(t / 50))
x <- rep(0, n)
for (i in 2:n) {
	x[i] <- x[i - 1] + mu[i - 1] * (t[i] - t[i - 1]) + sigma[i - 1] * sqrt(t[i]- t[i - 1]) * rnorm(1)
}
plot(t, x, type = "l")

# 4. Brownian motion Brownian bridge construction
# Given u<s<t, W(u), and W(t), W(s)|W(u),W(t)~N(((t-s)W(u)+(s-u)W(t))/(t-u), (s-u)(t-s)/(t-u)
t <- seq(0, 100, 0.5)
n <- length(t)
library(data.table) # use data.table to enable pass by reference
x <- data.table(rep(0, n))
x[n] <- 2
left <- 1
right <- n
f <- function(left, right) {
	if (left < right - 1) {
		mid <- floor((left + right) / 2)
		tmp <- ((t[right] - t[mid]) * x[left] + (t[mid] - t[left]) * x[right]) / (t[right] - t[left]) + sqrt((t[right] - t[mid]) * (t[mid] - t[left]) / (t[right] - t[left])) * rnorm(1)
		x[mid, V1:=tmp]
		f(left, mid)
		f(mid, right)
	}
}
f(left, right)
plot(t, x$V1, type = "l")

# 5. Geometric Browian motion
mu <- 0.05
sigma <- 0.2
t <- seq(0, 1, 1 / 252)
n <- length(t) - 1
x <- rep(0, n + 1)
for (i in 2:(n + 1)) {
	x[i] <- x[i - 1] + mu * (t[i] - t[i - 1]) + sigma * sqrt(t[i]- t[i - 1]) * rnorm(1)
}
S0 <- 100
S <- S0 * exp(x)
plot(t, S, type = "l")

# 6. CIR process
# dr=a(b-r)dt+sigma*sqrt(r)dW
t <- seq(0, 2, 1 / 252)
n <- length(t) - 1
b <- 1.2
a <- 1
sigma <- 0.2
r <- rep(0, n + 1)
r[1] <- 0.5
for (i in 2:(n + 1)) {
	r[i] <- r[i - 1] + a * (b - r[i - 1]) * (t[i] - t[i - 1]) + sigma * sqrt(max(r[i - 1], 0)) * sqrt(t[i] - t[i - 1]) * rnorm(1)
}
plot(t, r, type = "l")

# 7. Merton's jump diffusion model
# dS(t)/S(t-)=mu*dt+sigma*dW(t)+dJ(t), J(t)=sum_{j=1}^N(t)(Y(j)-1)
# N(t)~Poisson(lambda*t), Y(j)~LogNormal
# S(t)=S(0)e^{(mu-sigma&2/2)t+sigma*W(t)}\prod_{j=1}^N(t)Y(j)
mu <- 0.05
sigma <- 0.2
lmd <- 2
lmu <- 0.2
lsigma <- 0.4
t <- seq(0, 1, 1 / 252)
n <- length(t)
S0 <- 100
x <- rep(0, n)
for (i in 1:(n - 1)) {
	N <- rpois(1, lmd * (t[i + 1] - t[i]))
	M <- ifelse(N == 0, 0, sum(rnorm(N, lmu, lsigma)))
	x[i + 1] <- x[i] + (mu - sigma ^2 / 2) * (t[i + 1] - t[i]) + sigma * sqrt(t[i + 1] - t[i]) * rnorm(1) + M
}
S <- S0 * exp(x)
plot(t, S, type = "l")