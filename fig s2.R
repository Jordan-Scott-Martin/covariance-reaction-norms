setwd(...)

#generate values
b = rnorm(1e4, 0, 1)
x = seq(-1, 1, by = 0.1)
z = 1 + rep(b, each = length(x)) * exp(0.15 * x)
deltav = aggregate(z ~ rep(x, length(b)), FUN = var)
expv = var(b) * exp(0.3 * x) 

#observed variance
plot(deltav$z ~ x, type = "l", lwd = 3)
#expected variance = simulated variance
lines( expv ~ x, lwd = 3, lty = "dotted", col = "pink") 

#plot some individual responses
bi = c(0.1, 0.2, 0.3)
df = data.frame(
  z = 1 + rep(bi, each = length(x)) * exp(0.15 * x),
  x = rep(x, length(bi)),
  bi = rep(bi, each = length(x)) )

png("fig s2.png", width = 6, height = 6, units = "in", res = 300)
plot(z ~ x, type = "l", lwd = 3, col = "blue", data = df[df$bi==0.1,], ylim = c(1,1.5))
lines(z ~ x, type = "l", lwd = 3, col = "seagreen", data = df[df$bi==0.2,])
lines(z ~ x, type = "l", lwd = 3, col = "orange", data = df[df$bi==0.3,])
text(x = 0.9, y = 1.15, labels = expression(b[i] == 0.1), col = "blue", cex = 1.2)
text(x = 0.9, y = 1.26, labels = expression(b[i] == 0.2), col = "seagreen", cex = 1.2)
text(x = 0.9, y = 1.38, labels = expression(b[i] == 0.3), col = "orange", cex = 1.2)
text(x = -0.5, y = 1.48, labels = expression(z == 1 + b[i] * exp(0.15 * x)), col = "black", cex = 1.5)
dev.off()
