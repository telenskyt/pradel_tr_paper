# Function to prepare CES data for the model

source("functions/ces-data-capture_history.R")
source("functions/marray.R")

# Function to prepare CES data for the model
#
# input: 
# - ces - table (data.frame) with all individual captures. Following columns are assumed:
#		- exID - unique ID of the individual bird
#		- euring - species code
#		- site_ces - number code of the site
#		- year
#		- age - euring code of age (>= 4 means adult age)
#		
# - sp_id - species id
# - YRS - the years to consider
# - no_known_resid - if TRUE, it will "throw out" the known resid info (it will set num.captures.first.year = 1) - normally you don't want this, this is just for an experimentation
#
# output: a list with:
# 	- F, R: site-wise arrays of m-arrays
#	- k (number of temporal occasions)
#	- species_sites - list of CES sites where the species occurs
#	- visit: table of visits (row ~ site, column ~ year), equals 1 if visit was performed

prepare_data <- function(ces, sp_id, YRS, no_known_resid = FALSE)
{

	k <- n.occasions <- length(YRS)

	filter <- ces$euring == sp_id & ces$age >= 4 & ces$year %in% YRS # filter adults of given species
	sp_sites <- unique(sort(ces[filter,'site_ces']))
		# sites where species is, and it is at the same time a mapping to numbers 1..number of sites

	F <- array(data = NA, dim = c(length(sp_sites), k, k+1))
	R <- array(data = NA, dim = c(length(sp_sites), k-1, k))

	for (site_ind in 1:length(sp_sites)) {
		site <- sp_sites[site_ind]
		CH <- create_CH(ces, filter & ces$site_ces == site, YRS, NAs_in_missing_visits = FALSE)
		num.captures.first.year <- num_captures_first_year(ces, filter & ces$site_ces == site, YRS)
		if (no_known_resid)
			num.captures.first.year <- rep(1, length(num.captures.first.year))
		#known.resid <- as.integer(num.captures.first.year > 1)

		site_mar <- marray_multistage(CH, num.captures.first.year)
		F[site_ind,,] <- site_mar$F
		R[site_ind,,] <- site_mar$R
	}
	stopifnot(all(!is.na(F)))
	stopifnot(all(!is.na(R)))

	{ # build table of visits
		ces_yr_visits <- sqldf("
		select site_ces, year
		from ces
		group by site_ces, year
		") # global table, now restrict it to the sites relevant for our species
		ces_yr_visits$dummy <- TRUE
		stopifnot(!is.factor(ces_yr_visits$site_ces))
		visits1 <- acast(ces_yr_visits[ces_yr_visits$year %in% YRS & ces_yr_visits$site_ces %in% sp_sites,], site_ces ~ factor(year, levels = YRS), any, value.var = "dummy", drop = FALSE)
		stopifnot(all(as.integer(rownames(visits1)) == sp_sites))
		visit <- ifelse(visits1, 1, 0)
	}

	return(list(F = F, R = R, k = k, species_sites = sp_sites, visit = visit))
}

