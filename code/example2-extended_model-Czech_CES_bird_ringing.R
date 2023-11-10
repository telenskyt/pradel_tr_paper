#######
#
# Extension of the Pradel (1996) model for transients and multiple sites with different temporal coverage
#
# THIS IS THE MAIN SCRIPT TO RUN THE MODEL, START HERE

library(coda)

fast_model_run <- TRUE # very brief model run (few iterations) just for testing
run_parallel <- TRUE # disable only for debugging

source("../data/ces-data.R", chdir = TRUE) # load CES data

source("functions/ces-data-prepare.R") # functions to prepare CES data for the model

source("functions/nimble-wrapper.R") # functions to run Nimble

### prepare the data

YRS <- min(ces$year):max(ces$year) # year span
#sp_id <- 12510 # Euring species id, 12510 = Acrocephalus scirpaceus
sp_id <- 14640 # Euring species id, 14640 = Parus major

# process the individual captures and prepare the data in the required format:
data <- prepare_data(ces, sp_id, YRS)

### run the model

source("../model/model-pradel_with_transients-extended.nimble.R") # run the extended model, using Nimble

### post-hoc calculation of the derived parameters. Performed on the MCMC samples
### to ensure proper propagation of uncertainty of the primary model parameters to these derived quantities
### (for more explanation see e.g. Kéry, M., & Schaub, M. (2012). Bayesian Population Analysis using WinBUGS. Elsevier.)

source("functions/mcmc.R")

mcmcsum(out$mcmc) # summary of the MCMC chains

mc <- as.mcmc.list(out$mcmc)
mt <- as.matrix(mc, chains = TRUE, iters = TRUE)
n.occasions <- length(YRS)

# extract the vital rates calculated by the model

surv <- mt[,paste0('surv[',1:(n.occasions-1),']')]
sen <- mt[,paste0('sen[',1:(n.occasions-1),']')]

# calculate lambda = surv/sen

lambda <- matrix(nrow = nrow(mt), ncol = n.occasions - 1)
lambda <- surv/sen
colnames(lambda) <- paste0("lambda[", 1:ncol(lambda), "]")
mt <- cbind(mt, lambda)

# calculate recruitment

recr <- matrix(nrow = nrow(mt), ncol = n.occasions - 1)
recr <- lambda - surv
colnames(recr) <- paste0("recr", "[", 1:ncol(recr), "]")
mt <- cbind(mt, recr)

# calculate adult relative population index (equals 1 on the first occasion)

pop.idx <- matrix(nrow = nrow(mt), ncol = n.occasions)
pop.idx[,1] <- 1
for (i in 2:n.occasions) {
	pop.idx[,i] <- pop.idx[,i-1] * lambda[,i-1]
}
pop.idx.sum <- as.data.frame(cbind(Year = YRS, t(apply(pop.idx, 2, function (x) { c(Mean = mean(x), quantile(x, c(0.025, 0.975))) }))))

### As an example of a simple follow-up analysis with proper propagation of uncertainty, we calculate correlations
### of the demographic rates and density (population index). The follow-up analysis (here, correlation) is performed
### for each MCMC sample from the posterior distributions of the parameter estimates. This ensures proper propagation
### of uncertainty of the primary model parameters to the uncertainty of the resultant statistics (here, pearson correlation coefficient)
### (for more explanation see e.g. Kéry, M., & Schaub, M. (2012). Bayesian Population Analysis using WinBUGS. Elsevier.)

n.iter <- nrow(surv)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(surv[i,], recr[i,])$estimate
}
r.surv_recr <- c("r.surv_recr" = mean(r), quantile(r, c(0.025, 0.975)))

n.iter <- nrow(lambda)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(lambda[i,], surv[i,])$estimate
}
r.lambda_surv <- c("r.lambda_surv" = mean(r), quantile(r, c(0.025, 0.975)))

n.iter <- nrow(lambda)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(lambda[i,], recr[i,])$estimate
}
r.lambda_recr <- c("r.lambda_recr" = mean(r), quantile(r, c(0.025, 0.975)))

n.iter <- nrow(lambda)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(lambda[i,], pop.idx[i,1:(n.occasions-1)])$estimate
}
r.density_lambda <- c("r.density_lambda" = mean(r), quantile(r, c(0.025, 0.975)))

n.iter <- nrow(lambda)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(surv[i,], pop.idx[i,1:(n.occasions-1)])$estimate
}
r.density_surv <- c("r.density_surv" = mean(r), quantile(r, c(0.025, 0.975)))

n.iter <- nrow(lambda)
r <- c()
for (i in 1:n.iter) {
	r[i] <- cor.test(recr[i,], pop.idx[i,1:(n.occasions-1)])$estimate
}
r.density_recr <- c("r.density_recr" = mean(r), quantile(r, c(0.025, 0.975)))


# print the results
res <- rbind(r.surv_recr, r.lambda_surv, r.lambda_recr, r.density_lambda, r.density_surv, r.density_recr)
colnames(res)[1] <- 'Mean'
print(res)

##### Talon plot

library(talonplot)

big_font <- TRUE

talon_plot(
	recr_pts = colMeans(recr), surv_pts = colMeans(surv),
	recr_samples = recr, surv_samples = surv,
	recr_CI_low = apply(recr, 2, quantile, 0.025),
	recr_CI_high = apply(recr, 2, quantile, 0.975),
	surv_CI_low = apply(surv, 2, quantile, 0.025),
	surv_CI_high = apply(surv, 2, quantile, 0.975),
	CI_type = "ellipse",
	CI_transform = FALSE,
	main = "Parus major",
	big_font = big_font,
	plot.samples = TRUE,
	col.samples = "black"
)



