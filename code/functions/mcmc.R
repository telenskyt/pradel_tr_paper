
# summary of mcmc.list, with convergence stats
# input: mc - mcmc.list
mcmcsum <- function (mc, ...) 
{
	stopifnot(class(mc) == "mcmc.list")

	require(coda)

	mt <- as.matrix(mc)

	s1 <- summary(mc) 
	s2 <- cbind(s1$statistics[,c(1,2)], s1$quantiles[,c(1,5)], s1$statistics[,4, drop = FALSE])
	psrf <- tryCatch({
			gelman.diag(mc, transform=FALSE, autoburnin=FALSE, ...)
		},
		error = function (e) {
			psrf <- gelman.diag(mc, transform=FALSE, autoburnin=FALSE, multivariate = FALSE)
			psrf$mpsrf <- NA
			psrf
		})

	if (0) { # without psrf CI
		psrf.val <- psrf$psrf[,1,drop = FALSE]
		colnames(psrf.val) <- "psrf"
	} else {
		psrf.val <- psrf$psrf
		colnames(psrf.val) <- c("psrf", "psrf CI")
	}
	if (all(rownames(s2) == rownames(psrf.val))) { # this should normally be TRUE
		s4 <- cbind(s2, psrf.val)
	}
	else { # old way with merge - but it changes the order. It shouldn't be needed.
		s3 <- merge(s2, psrf.val, by = "row.names")
		row.names(s3) <- s3[,1]
		s4 <- s3[,-1]
	}

	list(summ = s4, mpsrf = psrf$mpsrf)
}