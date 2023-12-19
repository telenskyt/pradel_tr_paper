# Functions to take raw capture data and reshape them for the model with the two extensions:
#	- multiple sites with potentially different temporal coverage (Model extension 1 in the manuscript);
#	- confirmation of residency status within the first capture occasion (Model extension 2 in the manuscript).
# It basically creates a multi-site m-array.

source("functions/data-capture_history.R")
source("functions/marray.R")

# @description
# Prepare data for the model
# Prepare data for the extension of the Pradel (1996) model accounting for transients (Telensky et al. 2023, https://doi.org/10.1111/2041-210X.14262),
# in particular, for the extended variant, featuring:
#	- multiple sites with potentially different temporal coverage (Model extension 1 in the manuscript);
#	- confirmation of residency status within the first capture occasion (Model extension 2 in the manuscript).
#
# It basically creates a multi-site m-array.
#
# @param cap - table (data.frame) with individual captures (might be e.g. all captures of certain species) - one line per capture. 
#		Multiple captures allowed. Following columns are assumed:
#		- site - numeric code of the site
#		- exID - unique ID of the individual animal. If the same exID appears in multiple sites, it is treated as if it was a different individual.
#		- time - integer, the temporal occasion (e.g. a year); can start at any number; +1 means next occasion.
#
# @param visits - A data.frame with the columns `site` and `time` (same as in the `cap` parameter), 
#		one row for each site & time combination that waxtmation about missing visits, and these might differ between sites)

prepare_data <- function(cap, visits, known_resid_first_cap = NULL)
{
	occasions <- min(cap$time):max(cap$time) # occasions covering the capture data
	k <- n.occasions <- length(occasions)

	sp_sites <- unique(sort(cap[,'site']))
		# sites with some capture data; at the same time, it is a mapping to numbers 1..[number of sites]

	F <- array(data = NA, dim = c(length(sp_sites), k, k+1))
	R <- array(data = NA, dim = c(length(sp_sites), k-1, k))

	for (site_ind in 1:length(sp_sites)) {
		site <- sp_sites[site_ind]
		CH <- create_CH(cap[cap$site == site,], occasions = occasions, NAs_in_missing_visits = FALSE)
		if (!is.null(known_resid_first_cap)) {
			known.resid.first.cap <- known_resid_first_cap(cap[cap$site == site,])
			stopifnot(nrow(known.resid.first.cap) == nrow(CH)) # not only the count, but the rows must directly correspond!
		} else 
			known.resid.first.cap <- NULL

		site_mar <- marray_multistage(CH, known.resid.first.cap)
		F[site_ind,,] <- site_mar$F
		R[site_ind,,] <- site_mar$R
	}
	stopifnot(all(!is.na(F)))
	stopifnot(all(!is.na(R)))
	
	{ # build table of visits specific for the selection of captures which are in the `cap` data set
		visits$dummy <- TRUE
		stopifnot(!is.factor(visits$site))
		stopifnot(all(occasions %in% unique(visits$time)))
		visits1 <- acast(visits[visits$time %in% occasions & visits$site %in% sp_sites,], site ~ factor(time, levels = occasions), any, value.var = "dummy", drop = FALSE)
		stopifnot(all(as.integer(rownames(visits1)) == sp_sites))
		visit <- ifelse(visits1, 1, 0)
	}

	return(list(F = F, R = R, k = k, species_sites = sp_sites, visit = visit))
}

# @description
#  Returns a boolean vector that for each individual denotes whether it was recaptured within the occasion of the first capture, 
#  i.e. whether it was captured during multiple sub-occasions within the occasion of the first capture
#		(for more details see `Model extension 2: Confirmation of residency status within the first capture occasion` in the paper
#		 Telenský, T., Storch, D., Klvaňa, P., Reif, J. Extension of Pradel capture-recapture survival-recruitment model accounting for transients. 
#		 Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14262)
#
# @param cap - table (data.frame) with individual captures (might be e.g. all captures of certain species) - one line per capture. 
#		Multiple captures allowed. Following columns are assumed:
#		- exID - unique ID of the individual animal
#		- time - integer, the temporal occasion (e.g. a year); can start at any number; +1 means next occasion.
#		- column named as given by `visit.col` parameter - a unique identifier of a visit, i.e. a sub-occasion within the occasion (`time`)
#
# @param visit.col - name of the column in the `cap` table denoting the unique identifier of a visit (a sub-occasion)

recaptured_on_first_occasion <- function (cap, visit.col = "visit")
{
	CH0sum <- acast(cap, exID ~ time, function (x) length(unique(x)), value.var = visit.col)

	get.first <- function(x) min(which(x>0))
	num.captures.first.occasion <- apply(CH0sum, 1, function (x) x[get.first(x)])
	num.captures.first.occasion > 1
}

# @description Return a table with all combinations of site & time where visit was performed.
#
# @param cap - table (data.frame) with __ALL__ individual captures - one line per capture. 
#		Multiple captures allowed. Following columns are assumed:
#		- site - numeric code of the site
#		- time - integer, the temporal occasion (e.g. a year); can start at any number; +1 means next occasion.
#
# @details The function assumes that `cap` contains some capture record for every site & time combination where visit was performed.
#
# @return A data.frame with the columns `site` and `time` (same as in the `cap` parameter), 
#		one line for each site & time combination where the visit (capture occasion) was performed.
visits <- function(cap)
{
	sqldf("
	select site, time
	from cap
	group by site, time
	")
}

