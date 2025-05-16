
x = seq(-1, 1, by = 0.01)
vexp  = exp(x)
vquad  = x^2
vsp  = log(1 + exp(x))  # softplus
sdexp  = sqrt(exp(x))
sdquad  = x
sdsp  = sqrt(log(1 + exp(x)))  # softplus

col_exp  = "#ff395d"
col_quad  = "#541bff"
col_sp  = "#aa2aae"

png("link function comparison.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow=c(1,2))

#Variances
plot(vexp ~ x, ylim = c(0, 2.8), xlim = c(-1, 1), 
     ylab = expression(variance~(sigma^2)), 
     xlab = expression(linear~predictor~(eta)),
     type = "l", col = col_exp, lwd = 3)
lines(vquad ~ x, col = col_quad, lwd = 3)
lines(vsp ~ x, col = col_sp, lwd = 3)
legend("topleft", legend = c(
  expression(e^eta),
  expression(eta^2),
  expression(log(1 + e^eta))),
col = c(col_exp, col_quad, col_sp), lwd = 3, bty = "n")

#SDs
plot(sdexp ~ x, ylim = c(0, 2.8), xlim = c(-1, 1), 
     ylab = expression(standard~deviation~(sigma)), 
     xlab = expression(linear~predictor~(eta)),
     type = "l", col = col_exp, lwd = 3)
lines(sdquad ~ x, col = col_quad, lwd = 3)
lines(sdsp ~ x, col = col_sp, lwd = 3)
legend("topleft", legend = c(
  expression(sqrt(e^eta)),
  expression(eta),
  expression(sqrt(log(1 + e^eta)))),
col = c(col_exp, col_quad, col_sp), lwd = 3, bty = "n")

dev.off()
