######################################################
# Functions to transform capture histories to m-arrays
#
# Based on code from BPA 2012 (Kéry, M., & Schaub, M. (2012). Bayesian Population Analysis using WinBUGS. Elsevier.)


# marray()
# Create an m-array based on capture histories (CH)
# Code from BPA 2012, chapter 7.10 
#
# outputs m-array (lets denote as `mar'):
# k ........ number of occasions
# mar[i, t] ... i = 1..k-1; t = 1..k
# mar[i, t], for t = 1..k-1, is a number of animals first captured at occasion i next recaptured
#			at occasion t+1 (so recaptures within the same occasion are not counted here)
# mar[i, k] is number of animals first captured at occasion i and never recaptured again

marray <- function(CH) {
	nind <- nrow(CH)
	n.occasions <- ncol(CH)
	m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
	if (nind == 0) {
		out <- m.array[1:(n.occasions-1),2:(n.occasions+1)] # the same as below
		return(out)
	}
	# Calculate the number of released individuals at each time period
	for (t in 1:n.occasions) {
		m.array[t,1] <- sum(CH[,t])
	}
	for (i in 1:nind) {
		pos <- which(CH[i,]!=0)
		g <- length(pos)
		for (z in 1:(g-1)) {
			m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
		} #z
	} #i
	# Calculate the number of individuals that is never recaptured
	for (t in 1:n.occasions){
		m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
	}
	out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
	return(out)
}

# marray_multistage()
#
# Returns a pair of m-arrays F, R (F - for animals captured for the 1st time; R - for animals captured for 2nd and more-th time) 
#	for the purpose of transience
#
# Based on code from BPA 2012, chapter 7.10.3 (Kéry, M., & Schaub, M. (2012). Bayesian Population Analysis using WinBUGS. Elsevier.)
# 	  - there it is as an m-array with two age classes  - I adapted it for handling transience
# 	  analogy: juv -> F, then go to ad -> R; first captured as adults don't exist

# Returns: list with F, R, two m-arrays:
#	(k is number of temporal occasions)
# F[i,t] ... i = 1..k; t = 1..k+1
# 		F[i, t], for t = 1..k, is a number of animals first captured at occasion i next recaptured
#			at occasion t (t can be the same as i)
#		F[i, k+1] is number of animals first captured at occasion i and never recaptured again
# R[i, t] ... i = 1..k-1, t = 1..k
#		R[i, t], for t = 1..k-1, is a number of animals second captured at occasion i and next recaptured
#			at occasion t+1 (!)
#
# Note that m-arrays cannot handle different patterns in missing visits.

marray_multistage <- function (CH, num.captures.first.year)
{
	# input CH is without the transience trick
	n.occasions <- ncol(CH)

	# Transience trick (sometimes called extended capture histories): 
	# Behind the first capture, insert a known residency status (num.captures.first.year > 1), or NA if not known
	CH2 <- extend_CH(CH, num.captures.first.year) # output: CH2 - with the "transients trick" :)

	# Following data preparation part comes from the code in BPA chapter 7.10.3 
	CH.F <- CH2
	cap <- apply(CH.F, 1, sum)
	CH.F.R <- CH.F[cap >= 2,,drop=FALSE] # First captured CH recaptured at least once
	CH.F.N <- CH.F[cap <= 1,,drop=FALSE] # First captured CH never recaptured
	# Remove first capture
	first <- numeric()
	if (nrow(CH.F.R) > 0) {
		for (i in 1:nrow(CH.F.R)) {
			first[i] <- min(which(CH.F.R[i,]==1))
		}
	}
	CH.F.R1 <- CH.F.R
	if (nrow(CH.F.R) > 0) {
		for (i in 1:nrow(CH.F.R)) {
			CH.F.R1[i,first[i]] <- 0
		}
	}
	# Now this is a CH of captured again (i.e. residents); create m-array
	CH.R <- CH.F.R1
	CH.R.marray0 <- marray(CH.R)
	# And now my change!!!
	stopifnot(all(CH.R.marray0[1,] == 0))
	stopifnot(all(CH.R.marray0[,1] == 0))
	CH.R.marray <- CH.R.marray0[-1,-1]

	# Create CH matrix for first captured, ignoring subsequent recaptures
	second <- numeric()
	if (nrow(CH.F.R1) > 0) {
		for (i in 1:nrow(CH.F.R1)) {
			second[i] <- min(which(CH.F.R1[i,]==1))
		}
	}
	CH.F.R2 <- matrix(0, nrow = nrow(CH.F.R), ncol = ncol(CH.F.R))
	if (nrow(CH.F.R)) {
		for (i in 1:nrow(CH.F.R)){
			CH.F.R2[i,first[i]] <- 1
			CH.F.R2[i,second[i]] <- 1
		}
	}
	# Create m-array for these
	CH.F.R.marray <- marray(CH.F.R2)
	# The last column ought to show the number of first captured, not recaptured again
	# and should all be zeros, since all of them were recaptured
	CH.F.R.marray[,ncol(CH.F)] <- 0
	# Create the m-array for first captured (no juv., see the comment above) never recaptured and add it to the previous m-array
	CH.F.N.marray <- marray(CH.F.N)
	CH.F.marray <- CH.F.R.marray + CH.F.N.marray

	# data check:
	stopifnot(all(rowSums(CH.R.marray) == colSums(CH.F.marray[,1:(n.occasions-1)]) + 
						c(0, colSums(CH.R.marray[,1:(n.occasions-2)]))))

	return(list(F = CH.F.marray, R = CH.R.marray))
}
