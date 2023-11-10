##########################################################################################
#
# R-Nimble wrapper
# version 1: parallel setup, based on this example: https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html

run_nimble <- function (seed, code, bugs.data, bugs.constants = list(), dimensions = list(), inits, parameters, run_parallel = TRUE, ni, nt, nb, nc)
{

if (run_parallel) {

	# Parallel setup: based on https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html

	require(parallel)

	this_cluster <- makeCluster(nc)

	set.seed(234)

	clusterExport(this_cluster, varlist = c("ni", "nt", "nb", "nc"))
	output <- parLapply(cl = this_cluster, X = 1:nc,
							  fun = run_MCMC_allcode, 
							  code = code,
							  bugs.data = bugs.data,
							  bugs.constants = bugs.constants,
							  dimensions = dimensions,
							  inits = inits,
							  parameters = parameters, 
							  run_parallel = run_parallel,
							  ni = ni,
							  nt = nt,
							  nb = nb,
							  nc = nc) # similar as below
	# It's good practice to close the cluster when you're done with it.
	stopCluster(this_cluster)

	out <- list( 
		mcmc = as.mcmc.list(lapply(output, function (x) as.mcmc(x$runMCMC$samples))), # now I am using samplesAsCodaMCMC = TRUE, so the inner call as.mcmc() is not needed anymore, but its ok
		WAIC = mean(sapply(output, function (x) x$runMCMC$WAIC$WAIC)), # average over chains
		lppd = mean(sapply(output, function (x) x$runMCMC$WAIC$lppd)), # lppd: The log predictive density component of WAIC. (see ?waic)
		pWAIC = mean(sapply(output, function (x) x$runMCMC$WAIC$WAIC)), # pWAIC: as if "number of parameters" used for computation
		per_chain = list(
			WAIC = sapply(output, function (x) x$runMCMC$WAIC$WAIC),
			lppd = sapply(output, function (x) x$runMCMC$WAIC$lppd),
			pWAIC = sapply(output, function (x) x$runMCMC$WAIC$pWAIC)
		),
		seed = sapply(output, function (x) x$seed),
		ni = ni, nb = nb, nt = nt, nc = nc
	)
	return(out)
} else {
	warning("Non-parallel version doesn't yield results properly, use just for error debugging; it yields different output than the parallel version")
	output <- run_MCMC_allcode(seed = 234,
							  code = code,
							  bugs.data = bugs.data,
							  bugs.constants = bugs.constants,
							  dimensions = dimensions,
							  inits = inits,
							  parameters = parameters, 
							  run_parallel = FALSE,
							  ni = ni,
							  nt = nt,
							  nb = nb,
							  nc = nc) # similar as above

}
							  
}

run_MCMC_allcode <- function(seed, code, bugs.data, bugs.constants, dimensions = list(), inits, parameters, run_parallel = TRUE, ni, nt, nb, nc) 
{
	cat("Setting seed: \n")
	print(seed)
	cat("(end of seed)\n\n")

	set.seed(seed)

	require(nimble)

	myModel <- nimbleModel(code = code,
						  data = bugs.data,
						  constants = bugs.constants,
						  dimensions = dimensions,
						  inits = inits())

	CmyModel <- compileNimble(myModel)

	myMCMC <- buildMCMC(CmyModel, monitors = parameters, enableWAIC = TRUE)
	CmyMCMC <- compileNimble(myMCMC)

	if (run_parallel)
		nc <- 1
	results <- runMCMC(CmyMCMC, niter = ni*nt, nburnin = nb*nt, thin = nt, nchains = nc, setSeed = seed, WAIC = TRUE) # samplesAsCodaMCMC = TRUE - not needed any more, I am doing it myself

	gc()
	return(list(
		runMCMC = results,
		seed = seed
	))
}
