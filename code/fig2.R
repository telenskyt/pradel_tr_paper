# Figure 2
#
# This script generates a part of Figure 2 - a graph showing the increasing amount of recaptured individuals
# in a constant population.
#
# 

all_cap <- 1
rCt <- 0.7 # residency probability (within all captured animals)
p <- 0.4 # capture probability
sen <- 0.7 # seniority
T <- 5 # number of temporal occasions

# Calculate Equation 3 from the paper
xi <- c()  
xi[1] <- 1
for (i in 2:T) {
	xi[i] <- (1-sen) + sen*(1 - p)*xi[i-1]
}

# AC: allready captured before ("marked")
AC <- rCt * (1 - xi) 


par(cex = 1.3)
plot(1:T, AC, type = "l", ylim = c(0,1), xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n")
title(xlab = "time (occasions)", line = 2)
axis(1, pos = 0)
axis(2, at = c(0, rCt, 1), labels = c(0, "", 1), las = 1, pos = 1)
lines(c(1,T), c(rCt, rCt))
lines(c(1,T), c(1, 1))
polygon(c(1:T,T,1), c(AC,1,1), density = 10, angle = 45)
polygon(c(1:T,T,1), c(AC,0,0), density = 10, angle = -45)

savePlot("res_tr_graf", "wmf")
savePlot("res_tr_graf", "emf")

#abline(h = rCt) # this plots beyond the vertical axes, unfortunatelly
#abline(h = 1)


	