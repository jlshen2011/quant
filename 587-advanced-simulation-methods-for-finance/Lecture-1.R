# Simulation Methods in Financial Engineering
# R ode for Lecture 1
# Simulation from univariate distributions
# Jieli Shen


# 1. Simulation from a distribution

# Simulation from an exponential distribution
lambda <- 1
u <- runif(1000)
x <- -1 / lambda * log(u)
hist(x, freq = F)

# Simulation from a standard normal distribution
# Method 1: central limit theorem
u <- array(runif(12000), dim = c(1000, 12))
x <- apply(u, 1, sum) - 6
hist(x, freq = F)
# Method 2: Box-Muller transform
u1 <- runif(500)
u2 <- runif(500)
x <- c(sqrt(-2 * log(u2)) * cos(2 * pi * u1), sqrt(-2 * log(u2)) * sin(2 * pi * u1))
hist(x, freq = F)

# Simulation from a binomial distributions
p <- 0.3
n <- 10
u <- array(runif(1000 * n), dim = c(1000, n))
x <- apply(u, 1, function(x) sum(x < p))
hist(x, freq = F)

# Simulation from a Poisson distribution
# Method 1: use its relation to the exponential distribution
# x~Poisson(lambda) is the number of events in [0,1] when the time between consecutive events are i.i.d exp(lambda); then x is the integer k such that t1+...tk<=1 but t1+...+t(k+1)>1
lambda <- 5
x <- numeric(1000)
for (i in 1:1000) {
	s <- 0
	k <- 0
	while (s <= 1) {
		u <- runif(1)
		z <- -log(u) / lambda
		s <- s + z
		k <- k + 1
	}
	x[i] <- k - 1	
}
hist(x, freq = F)
# Method 2: inverse transform
lambda <- 5
x <- numeric(1000)
for (i in 1:1000) {
	k <- 0
	f <- exp(-lambda)
	F <- f
	u <- runif(1)
	while (F <= u){
		k <- k + 1
		f <- f * lambda / k
		F <- F + f
	}
	x[i] <- k
}
hist(x, freq = F)

# Simulation from a noncentral chi distribution
# chi^2_k(lambda)=(z1+lambda)^2+z2^2+...+zk^2
lambda <- 1
k <- 5
z <- array(rnorm(1000 * k), dim = c(1000, k))
z[, 1] <- z[, 1] + lambda
x <- apply(z, 1, function(x) sum(x ^ 2))
hist(x, freq = F)

# Simulation from a double exponential distribution
# f(x)=exp(-|x-mu|/theta)/(2*theta)
# Method 1: inveser transform
mu <- 0
theta <- 1
u <- runif(1000)
x <- mu - theta * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
hist(x, freq = F)
# Method 2: use its relation to the exponential distribution
u1 <- runif(1000)
u2 <- runif(1000)
z <- -theta * log(u1)
x <- ifelse(u2 < 0.5, mu + z, mu - z) 
hist(x, freq = F)