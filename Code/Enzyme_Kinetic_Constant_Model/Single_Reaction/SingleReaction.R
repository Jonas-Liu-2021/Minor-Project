library(deSolve)
library(ggplot2)

# Define the differential equation (Michaelis-Menten equation)
dSdt <- function(t, S, parms) {
  V_max <- parms[1]
  K_m <- parms[2]
  v <- - (V_max * S) / (K_m + S)
  return(list(v))
}

# Set the initial parameters (Vmax, Km, initial concentration of substrate)
V_max <- 0.02101
K_m <- 1.92
S0 <- 0.8
times <- seq(0, 700, length.out = 1000)

# Parameter adjustment
delta <- 0.05

# Solve the ODE for the different parameter adjustments
sol1 <- ode(y = S0, times = times, func = dSdt, parms = c(V_max, K_m))
sol2 <- ode(y = S0, times = times, func = dSdt, parms = c(V_max + delta * V_max, K_m))
sol3 <- ode(y = S0, times = times, func = dSdt, parms = c(V_max - delta * V_max, K_m))
sol4 <- ode(y = S0, times = times, func = dSdt, parms = c(V_max, K_m + delta * K_m))
sol5 <- ode(y = S0, times = times, func = dSdt, parms = c(V_max, K_m - delta * K_m))

# Calculate the absolute difference between the curves
diff1 <- abs(sol1[, 2] - sol2[, 2])
diff2 <- abs(sol1[, 2] - sol3[, 2])
diff3 <- abs(sol1[, 2] - sol4[, 2])
diff4 <- abs(sol1[, 2] - sol5[, 2])

# Plot the concentration-time curves
windows()
par(mfrow=c(1,2))
matplot(sol1[,1],sol1[,2],type = "l",col = 1, ylim = c(0,S0), xlab = "Time (minutes)",ylab = "Substrate Concentration (mM)", main="Concentration-Time Curves")
lines(sol2[,1],sol2[,2],type = "l",col=2)
lines(sol3[,1],sol3[,2],type = "l",col=3)
lines(sol4[,1],sol4[,2],type = "l",col=4)
lines(sol5[,1],sol5[,2],type = "l",col=5)
legend("topright",c("Unchanged", "Increase Vmax", "Decrease Vmax", "Increase Km", "Decrease Km"),lty = 1,col = 1:5,xjust = 0, bty="n",x.intersp=1)

# Plot the difference curves
matplot(sol1[,1],diff2, type = "l",xlab = "Time (minutes)",ylab = "Difference of Substrate Concentration (mM)", main="Difference Curves")
lines(sol1[,1],diff1, col=2)
lines(sol1[,1],diff3, col=3)
lines(sol1[,1],diff4, col=4)
legend("bottom",c("Increase Vmax", "Decrease Vmax", "Increase Km", "Decrease Km"),lty = 1,col = 1:4,xjust = 0, bty="n",x.intersp=1)

# Add points and coordinate labels to the maximum difference of each curve
max_index1 <- which.max(diff1)
max_diff1 <- round(max(diff1), 4)
points(times[max_index1], max_diff1, col = "red", pch = 19)
text(times[max_index1], max_diff1, labels = paste("(", round(times[max_index1], 2), ",", max_diff1, ")", sep = ""), pos = 4, col = "blue")

max_index2 <- which.max(diff2)
max_diff2 <- round(max(diff2), 4)
points(times[max_index2], max_diff2, col = "green", pch = 19)
text(times[max_index2], max_diff2, labels = paste("(", round(times[max_index2], 2), ",", max_diff2, ")", sep = ""), pos = 4, col = "blue")

max_index3 <- which.max(diff3)
max_diff3 <- round(max(diff3), 4)
points(times[max_index3], max_diff3, col = "blue", pch = 19)
text(times[max_index3], max_diff3, labels = paste("(", round(times[max_index3], 2), ",", max_diff3, ")", sep = ""), pos = 4, col = "blue")

max_index4 <- which.max(diff4)
max_diff4 <- round(max(diff4), 4)
points(times[max_index4], max_diff4, col = "purple", pch = 19)
text(times[max_index4], max_diff4, labels = paste("(", round(times[max_index4], 2), ",", max_diff4, ")", sep = ""), pos = 3, col = "blue")

