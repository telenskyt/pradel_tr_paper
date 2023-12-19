##########################################################################################
# Extension of the Pradel (1996) model accounting for transients (Telensky et al.)
# This is the simplest, basic model variant, estimating demographic parameters (survival, seniority, population growth, recruitment) over time
#
# Parametrization (this can be changed easily):
# - demographic parameters: estimate for each time occasion
# - capture probability: constant
# - residence probability: temporal slope (linear on the logit scale)
#
# Model written in BUGS, fit in Nimble.
#
# Script takes input from these global variables:
#	- data - result of call to prepare_data() - a list containing k, F, R (please see the detailed descriptions of these in the commentary of the model below, under "Input data")
#	- fast_model_run - FALSE for full model run, TRUE for just very fast test
# 	- run_parallel - TRUE/FALSE
#
# References: Telenský, T., Storch, D., Klvaňa, P., Reif, J. Extension of Pradel capture-recapture survival-recruitment model accounting for transients. Methods in Ecology and Evolution.

require(nimble)
require(coda)

source("../code/functions/nimble-wrapper.R") # functions to run Nimble
source("../code/functions/tools.R")

code <- nimbleCode({
#######################
#
# Model parameters:
#
# surv[t] .............. survival probability, probability that a resident animal present on occasion t will be present on occasion t + 1
#                        (meaning being available at the same site; survives "for the site", aka apparent survival)
# sen[t-1] ............. seniority probability, probability that a resident animal present on occasion t was already present on occasion t − 1 
# p[t] ................. capture probability, probability that resident animal alive and present on occasion t is captured on occasion t
# pr.captured_twice .... probability that resident animal is captured second time on the occasion of the first capture
# residency[t] ......... residency probability (within all captured); probability that animal captured on occasion t is a resident
# residency.first[t] ... residency probability (within first captured individuals); probability that animal captured first on occasion t is a resident		
#
# Input data (please find the detailed definitions in the manuscript, Table 1):#
# k .................... number of occasions
# F, R ................. two m-arrays:
#
# 	F[i,t] ... i = 1..k; t = 1..k+1
# 			F[i, t], for t = 1..k, t >= i, is a number of animals first captured at occasion i, next recaptured at occasion t
#			F[i, k+1] is number of animals first captured at occasion i and never recaptured again
# 	R[i, t] ... i = 1..k-1, t = 1..k
#			R[i, t], for t = 1..k-1, is a number of animals captured (not for the first time; i.e. confirmed residents) 
#				at occasion i and next recaptured at occasion t+1
#			R[i, k] is number of animals captured (not for the first time; i.e. confirmed residents) 
#				at occasion i and never recaptured again
# 
# F.i[t] ......... number of animals first captured at occasion t (row sums of F)
#
# time[i] .............. standardized number of the temporal occasion, used to calculate temporal slope
#
#######################


#### Priors

p.mean ~ dunif(0, 1) 

residency.mean ~ dunif(0, 1) # probability scale. Per habitat, 1 - dry, 2 - wet (reedbed, wet scrubland)
residency.time.slope ~ dnorm(0, 0.01)

for (i in 1:(k-1)) {
	surv[i] ~ dunif(0, 1)
	sen[i] ~ dunif(0, 1)
}

# Particularities of the model parameters 

for (i in 1:k) {
	logit(residency[i]) <- logit(residency.mean) + residency.time.slope * time[i] # temporal slope of residency
	p[i] <- p.mean
}

# Likelihood

# 1) Likelihood for row sums of F - the probability of being first captured on occasion t
p.resid.never.captured.before[1] <- 1
lambda_since_beg[1] <- 1
for (t in 2:k) {
	p.resid.never.captured.before[t] <- 1 - sen[t-1] + sen[t-1] * (1 - p[t-1]) * p.resid.never.captured.before[t-1] 
		# equation 3 in the manuscript
		# p.resid.never.captured.before - probability that resident animal present on occasion t was never captured on previous occasions; 
		#		or, equivalently, probability that resident animal captured on occasion t was captured for the first time

	lambda_since_beg[t] <- lambda_since_beg[t-1]* surv[t-1] / sen[t-1]
		# product of lambdas from temporal occasion 1 to t-1; that's in fact adult population index on occasion t
}
for (t in 1:k) {
	residency.first[t] <- p.resid.never.captured.before[t] * residency[t] / (p.resid.never.captured.before[t] * residency[t] + (1 - residency[t])) 
		# equation 4 in the manuscript
	
	f.i[t] <- p.resid.never.captured.before[t] * p[t] * lambda_since_beg[t] / residency.first[t]
}
f.i.normalized[1:k] <- f.i[1:k]/sum(f.i[1:k])
F.i[1:k] ~ dmulti(f.i.normalized[1:k], ex_caught) # division by the sum is basically the conditioning by being captured in any of the occasions

# 2) Likelihood for the individial F[i,t] - the probability of being captured second time given the first capture (thus becoming confirmed resident)
for (i in 1:k) {	 
	# compute f.i.t[i,t] - conditional on f.i[i] (F.i[i]) - probability that animal first caught at occasion i is next recaptured at occasion t
	for (t in 1:i) { # below diagonal, and the diagonal itself
		f.i.t[i,t] <- 0
	}
	for (t in (i+1):min(i+1,k)) { # need to do it separately from the following loop, because prod() cannot be empty..
		f.i.t[i,t] <- residency.first[i] * surv[i] * p[t]
	}
	for (t in (i+2):k) {
		f.i.t[i,t] <- residency.first[i] * surv[i] * prod(surv[(i+1):(t-1)]*(1 - p[(i+1):(t-1)])) * p[t]
	}
	f.i.t[i,k+1] <- 1 - sum(f.i.t[i,1:k])
	
	# Finally the probabilities for F
	F[i,1:(k+1)] ~ dmulti(f.i.t[i,], F.i[i])
}
# 3) Likelihood for R - probability of next recaptures, given the animal is a confirmed resident 
for (i in 1:(k-1)) {
	for (t in 1:(i-1)) { # below diagonal
		r.i.t[i,t] <- 0
	}
	r.i.t[i,i] <- surv[i] * p[i+1] # need to do it separately from the following, because prod() cannot be empty..
	for (t in (i+1):(k-1)) {
		# probability that animal captured (not for the first time; i.e. confirmed residents) 
		# at occasion i is next recaptured at occasion t+1
		r.i.t[i,t] <- surv[i] * prod(surv[(i+1):t]*(1 - p[(i+1):t])) * p[t+1]
	}
	r.i.t[i,k] <- 1 - sum(r.i.t[i,1:(k-1)])
	
	# Probabilities for R
	R[i,1:k] ~ dmulti(r.i.t[i,1:k], R.i[i])
}


#####
# Derived demographic parameters are calculated ex-post in the analysis script
# (it can be done incl. traceplots and psrf stats, it saves the size of the saved model)

})

# requires input: `data` (see comments above for description)

bugs.data <- list(
	F = data$F, 
	R = data$R,
	F.i = rowSums(data$F), 
	R.i = rowSums(data$R), 
	ex_caught = sum(data$F) # total number of animals caught
)
bugs.constants <- list(
	k = data$k,
	time = scale(1:data$k)[,1] # standardized index of temporal occasion (for use with temporal slope)
)
k <- bugs.constants$k
dimensions <- list( # Nimble sometimes requires dimensions of arrays to be specified
	f.i.t = c(k, k+1), 
	r.i.t = c(k-1, k),
	surv = k-1,
	sen = k-1
)
# Initial values
inits <- (function() {  # closure needs to be created, so that the variable k is available
			k <- k
			function(){
				list(p.mean = runif(1, 0.2, 0.25), residency.time.slope = runif(1, -0.5, 0.5), surv = runif(k-1, 0, 1), sen = c(runif(k-1, 0, 1)))
			}
		 })()

# Parameters monitored
parameters <- c("p.mean", "residency", "residency.mean", "residency.time.slope", "surv", "sen")

# MCMC settings
if (fast_model_run) { # brief model run just for testing
ni <- 100
nt <- 1
nb <- 50
nc <- 3
adapt <- 0
} else { # final run
ni <- 5000
nt <- 40
nb <- 2000
nc <- 3
adapt <- 5000
}

cat("ni = ", ni, ", nt = ", nt, ", nb = ", nb, ", nc = ", nc, "\n")

# Run the model in Nimble

mstart(absolute_time = TRUE)

out <- run_nimble(seed, code, bugs.data, bugs.constants, dimensions, inits, parameters, run_parallel = run_parallel, ni, nt, nb, nc)
							  
mstop(absolute_time = TRUE)

