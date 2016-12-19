# Simulation Methods in Financial Engineering
# R ode for Lecture 5
# Jieli Shen


# 1. Simulation from multivariate Gaussian distribution
mu <- c(-2, 1)
Sigma <- array(c(1, 0.5, 0.5, 2), dim = c(2, 2))
L.prime <- chol(Sigma)
L <- t(L.prime)
z <- array(rnorm(2 * 1000), dim = c(1000, 2))
x <- t(apply(z, 1, function(x) mu + L %*% x))
plot(x[, 1], x[, 2], xlab = expression(x[1]), ylab = expression(x[2]))


# 2. Algorithm for Gaussian copula
#	Simulate z_1,...,z_p ~ N(0,1)
#	Set z^*=Az ~ N(0,Sigma) where Sigma=AA^t
#	For i=1,...,p, set u_i=Phi(z_i^*/sigma_i), then u_i~U(0,1); set x_i=F_i^{-1}(u_i), 
#	then x_i~F_i
x <- array(0, dim = c(1000, 2))
for (i in 1:1000) {
	z <- rnorm(2)
	z.star <- L %*% z
	u <- c(pnorm(z.star[1] / sqrt(Sigma[1, 1])), pnorm(z.star[2] / sqrt(Sigma[2, 2])))
	x[i, ] <- c(qchisq(u[1], df = 10), qchisq(u[2], df = 3))
}
plot(x[, 1], x[, 2], xlab = expression(x[1]), ylab = expression(x[2]))