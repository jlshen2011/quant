# Simulation Methods in Financial Engineering
# R ode for Lecture 3
# Jieli Shen


# 1. Gibbs sampler


# 2. Metropolis-Hastings
# Algorithm:
# 	Pick initial values x(0)
# 	At each iteration k, choose a new x(k+1) as follows:
# 		Sample z~g(|x(k))
#		Accept z as x(k+1) with probability alpha=min(1, f(z)g(x(k)|z)/f(x(k))g(z|x(k))))
#		Otherwise set x(k+1)=x(k)
# Example:
#	Simulation from f(x)=f_t5(x)(1-sin(20x))/c for |x|<3
x <- numeric(10000)
x[1] <- 0
for (i in 2:10000) {
	z <- rt(1, df = 5)
	alpha <- min(1, (1 - sin(20 * z)) / (1 - sin(20 * x[i - 1])) * (abs(z) < 3))
	if (runif(1) < alpha) {
		x[i] <- z
	}
	else{
		x[i] <- x[i - 1]
	}
}
hist(x[5001:10000], freq = F, breaks = 100)