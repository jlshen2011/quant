# Simulation Methods in Financial Engineering
# R ode for Lecture 2
# Rejection sampling and importance sampling
# Jieli Shen


# 1. Rejection sampling
# Simulate N(0,1) from a double exponential distribution using rejection sampling
# N(0,1): f(x)=exp(-x^2/2)/sqrt(2*pi)
# DE: g(x)=exp(-abs(x))/2
# f(x)/g(x)=sqrt(2/pi)exp(-x^2/2+abs(x))<sqrt(2/pi)*exp(1/2)=1.32
f <- function(x) exp(-x ^ 2 / 2) / sqrt(2 * pi)
g <- function(x) exp(-abs(x)) / 2
const <- 1.3155
G.inv <- function(u) ifelse(u < 1 / 2, log(2 * u), -log(2 * (1 - u)))
i <- 0
x <- NULL
while (i <= 1000) {
	u1 <- runif(1)
	z <- G.inv(u1)
	u2 <- runif(1)
	if (u2 < f(z) / (const * g(z))) {
		x <- c(x, z)
		i = i + 1
	}
}
hist(x, freq = F)


# 2. Importance sampling
# Draw a large number of z~g(x)
# Sample a value x from {z_1,...,z_n} according to weights w_i=f(z_i)/g(z_i)
u <- runif(10000)
z <- G.inv(u)
w <- f(z) / g(z)
w <- cumsum(w) / sum(w)
x <- numeric(1000)
for (i in 1:1000) {
	u <- runif(1)
	x[i] <- z[sum(w < u) + 1]
}
hist(x, freq = F)