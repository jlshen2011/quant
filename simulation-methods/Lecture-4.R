# Simulation Methods in Financial Engineering
# R ode for Lecture 4
# Bootstrap methods
# Jieli Shen


ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
mydata <- data.frame(weight, group)

# Linear regression
lm.D9 <- lm(weight ~ group, data = mydata)
summary(lm.D9)
confint(lm.D9, "groupTrt")

# Model-based boostrap
boot.meth1 <- function(data = mydata, indices){
	data <- data[indices, ]
	mod <- lm(weight ~ group, data = data)
	coefficients(mod)
}

# Residual boostrap
boot.meth2 <- function(data = mydata, indices, fit = lm.D9){
	weight.boot <- fitted(lm.D9) + residuals(lm.D9)[indices]
	data.star <- data
	data.star[,1] <- weight.boot
	mod <- lm(weight ~ group, data = data.star)
	coefficients(mod)
}

# Use R's boot package
library(boot)
out.boot.meth1 <- boot(mydata, boot.meth1, 5000)
out.boot.meth1
boot.ci(out.boot.meth1, index = 2, type = c("norm", "perc"))
out.boot.meth2 <- boot(mydata, boot.meth2, 5000)
out.boot.meth2
boot.ci(out.boot.meth2, index = 2, type = c("norm", "perc"))
plot(out.boot.meth1, index = 2)
plot(out.boot.meth2, index = 2)