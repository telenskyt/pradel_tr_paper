##########################################################################################
# Extension of the Pradel (1996) model accounting for transients (Telensky et al.)
# This is an extended model variant, featuring:
#	- multiple sites with potentially different temporal coverage (Model extension 1 in the manuscript);
#	- confirmation of residency status within the first capture occasion (Model extension 2 in the manuscript).
# 
# Parametrization (this can be changed easily):
# - demographic parameters: estimate for each time occasion
# - capture probability and omega: constant over time, random site effects
# - residence probability: temporal slope (linear on the logit scale)
#
# Model written in BUGS, fit in Nimble.
#
# Script takes input from these global variables:
#	- data - result of call to prepare_data() - a list containing k, F, R, visit (please see the detailed descriptions of these in the commentary of the model below, under "Input data")
#	- fast_model_run - FALSE for full model run, TRUE for just very fast test
# 	- run_parallel - TRUE/FALSE
#
# References: Telenský, T., Storch, D., Klvaňa, P., Reif, J. Extension of Pradel capture-recapture survival-recruitment model accounting for transients. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14262

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
# 	F[site,i,t] ... i = 1..k; t = 1..k+1
# 			F[site, i, t], for t = 1..k, t >= i, is a number of animals first captured at occasion i, next recaptured at occasion t
#			F[site, i, k+1] is number of animals first captured at occasion i and never recaptured again
# 	R[site, i, t] ... i = 1..k-1, t = 1..k
#			R[site, i, t], for t = 1..k-1, is a number of animals captured (not for the first time; i.e. confirmed residents) 
#				at occasion i and next recaptured at occasion t+1
#			R[site, i, k] is number of animals captured (not for the first time; i.e. confirmed residents) 
#				at occasion i and never recaptured again
# 
# F.i[site, t] ......... number of animals first captured at occasion t (row sums of F)
#
# visit[site, i] ....... 0 - there was no visit at site on occasion i; 1 - there was a visit
#
# time[i] .............. standardized number of the temporal occasion, used to calculate temporal slope
#
#######################


#### Priors

p.mean ~ dunif(0, 1) 
p.mean_value <- mean(ilogit(logit(p.mean) + p.site.reff[1:sites])) # mean of the actual p[site], on the probability scale
p.site.sigma ~ dunif(0, 10)

pr.captured_twice.mean ~ dunif(0, 1) # jakoby intercept
pr.captured_twice.mean_value <- mean(pr.captured_twice[1:sites]) # mean of the actual pr.captured_twice[site], on the probability scale
pr.captured_twice.site.sigma ~ dunif(0, 10)

residency.mean ~ dunif(0, 1) # probability scale. Per habitat, 1 - dry, 2 - wet (reedbed, wet scrubland)
residency.time.slope ~ dnorm(0, 0.01)

for (i in 1:(k-1)) {
	surv[i] ~ dunif(0, 1)
	sen[i] ~ dunif(0, 1)
}

# Particularities of the model parameters 

for (site in 1:sites) {
	p.site.reff[site] ~ dnorm(0, 1/p.site.sigma^2) # site random effect for p
	pr.captured_twice.site.reff[site] ~ dnorm(0, 1/pr.captured_twice.site.sigma^2) # site random effect for pr.captured_twice
	logit(pr.captured_twice[site]) <- logit(pr.captured_twice.mean) + pr.captured_twice.site.reff[site]
}

for (i in 1:k) {
	logit(residency.time[i]) <- logit(residency.mean) + residency.time.slope * time[i] # temporal slope of residency
	for (site in 1:sites) {
		logit(residency[site, i]) <- logit(residency.time[i])
		p[site, i] <- visit[site, i] * ilogit(logit(p.mean) + p.site.reff[site]) # handling visits + random site effect
	}
}

# Likelihood
for (site in 1:sites) {
	# 1) Likelihood for row sums of F - the probability of being first captured on occasion t
	p.resid.never.captured.before[site, 1] <- 1
	lambda_since_beg[site, 1] <- 1
	for (t in 2:k) {
		p.resid.never.captured.before[site, t] <- 1 - sen[t-1] + sen[t-1] * (1 - p[site, t-1]) * p.resid.never.captured.before[site, t-1] 
			# equation 3 in the manuscript
			# p.resid.never.captured.before - probability that resident animal present on occasion t was never captured on previous occasions; 
			#		or, equivalently, probability that resident animal captured on occasion t was captured for the first time

		lambda_since_beg[site, t] <- lambda_since_beg[site, t-1]* surv[t-1] / sen[t-1]
			# product of lambdas from temporal occasion 1 to t-1; that's in fact adult population index on occasion t
	}
	for (t in 1:k) {
		residency.first[site, t] <- p.resid.never.captured.before[site, t] * residency[site, t] / (p.resid.never.captured.before[site, t] * residency[site, t] + (1 - residency[site, t])) 
			# equation 4 in the manuscript
		
		f.i[site, t] <- p.resid.never.captured.before[site, t] * p[site, t] * lambda_since_beg[site, t] / residency.first[site, t]
	}
	f.i.normalized[site,1:k] <- f.i[site,1:k]/sum(f.i[site,1:k])
	F.i[site,1:k] ~ dmulti(f.i.normalized[site,1:k], ex_caught[site]) # division by the sum is basically the conditioning by being captured in any of the occasions

	# 2) Likelihood for the individial F[site,i,t] - the probability of being captured second time given the first capture (thus becoming confirmed resident)
	for (i in 1:k) {	 
		# compute f.i.t[site,i,t] - conditional on f.i[site, i] (F.i[site, i]) - probability that animal first caught at occasion i is next recaptured at occasion t
		for (t in 1:(i-1)) { # below diagonal
			f.i.t[site,i,t] <- 0
		}
		f.i.t[site,i,i] <- residency.first[site,i] * pr.captured_twice[site]
		for (t in (i+1):min(i+1,k)) { # need to do it separately from the following loop, because prod() cannot be empty..
			f.i.t[site,i,t] <- residency.first[site,i] * (1 - pr.captured_twice[site]) * surv[i] * p[site,t]
		}
		for (t in (i+2):k) {
			f.i.t[site,i,t] <- residency.first[site,i] * (1 - pr.captured_twice[site]) * surv[i] * prod(surv[(i+1):(t-1)]*(1 - p[site,(i+1):(t-1)])) * p[site,t]
		}
		f.i.t[site,i,k+1] <- 1 - sum(f.i.t[site,i,1:k])
		
		# Finally the probabilities for F
		F[site,i,1:(k+1)] ~ dmulti(f.i.t[site,i,], F.i[site,i])
	}
	# 3) Likelihood for R - probability of next recaptures, given the animal is a confirmed resident 
	for (i in 1:(k-1)) {
		for (t in 1:(i-1)) { # below diagonal
			r.i.t[site,i,t] <- 0
		}
		r.i.t[site,i,i] <- surv[i] * p[site,i+1] # need to do it separately from the following, because prod() cannot be empty..
		for (t in (i+1):(k-1)) {
			# probability that animal captured (not for the first time; i.e. confirmed residents) 
			# at occasion i is next recaptured at occasion t+1
			r.i.t[site,i,t] <- surv[i] * prod(surv[(i+1):t]*(1 - p[site,(i+1):t])) * p[site,t+1]
		}
		r.i.t[site,i,k] <- 1 - sum(r.i.t[site,i,1:(k-1)])
		
		# Probabilities for R
		R[site,i,1:k] ~ dmulti(r.i.t[site,i,1:k], R.i[site,i])
	}
}

#####
# Derived demographic parameters are calculated ex-post in the analysis script
# (it can be done incl. traceplots and psrf stats, it saves the size of the saved model)

})

# input: data

bugs.data <- list(
	F = data$F, 
	R = data$R,
	F.i = apply(data$F, c(1,2), sum), 
	R.i = apply(data$R, c(1,2), sum), 
	ex_caught = apply(data$F, 1, sum), # total number of animals caught for each site
	visit = data$visit
)
bugs.constants <- list(
	k = data$k,
	sites = length(data$species_sites),
	time = scale(1:data$k)[,1] # standardized index of temporal occasion (for use with temporal slope)
)
k <- bugs.constants$k
sites <- bugs.constants$sites	
dimensions <- list( # Nimble sometimes requires dimensions of arrays to be specified
	f.i.t = c(sites, k, k+1), 
	r.i.t = c(sites, k-1, k),
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
parameters <- c("p.mean", "p.mean_value", "p.site.sigma", "residency.time", "residency.mean", "residency.time.slope", "pr.captured_twice.mean", "pr.captured_twice.mean_value", "pr.captured_twice.site.sigma", "surv", "sen")

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