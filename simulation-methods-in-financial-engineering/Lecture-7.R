# Simulation Methods in Financial Engineering
# R ode for Lecture 7
# Variance reduction techinques
# Jieli Shen


# 1. Control variates
# Estimate EX using mean(X)-cov(X,Y)/var(X)*(mean(X)-EX)
# Estimate Ef(x), f(x)=1/(1+x), x in [0,1] using g(x)=1+x as control variates
f1.mean <- f2.mean <- numeric(1000)
for (i in 1:1000) {
	u <- runif(1000)
	f1.mean[i] <- mean(1 / (1 + u))
	f2.mean[i] <- mean(1 / (1 + u)) - cov(1 / (1 + u), 1 + u) / var(1 + u) * (mean(1 + u) - 1.5)
}
c(mean(f1.mean), mean(f2.mean))
c(var(f1.mean), var(f2.mean))


# 2. Antithetic variates
lambda <- 1
x.mean <- y.mean <- z.mean <- numeric(1000)
for (i in 1:1000) {
	u1 <- runif(1000, 0, 1)
	x <- -1 / lambda * log(u1)
	x.mean[i] <- mean(x)
	u2 <- 1 - u1
	y <- -1 / lambda * log(u2)
	y.mean[i] <- mean(c(x, y))
	u3 <- runif(2000, 0, 1)
	z <- -1 / lambda * log(u3)
	z.mean[i] <- mean(z)
}
c(var(x.mean), var(y.mean), var(z.mean))
# [1] 0.0010544000 0.0001725821 0.0005181840
# If x and y are negatively correlated, this is even better than using 2n i.i.d. samples


# 3. Stratified sampling
# Sampling each strata Nj samples proportional to pj for all j
x1 <- x2 <- numeric(1000)
for(i in 1:1000) {
	x1[i] <- mean(runif(100))
	for(j in 1:10) {
		x2[i] = x2[i] + 0.1 * mean(runif(10, (j - 1) * 0.1, j * 0.1))
	}
}
c(mean(x1), mean(x2))
c(var(x1), var(x2))