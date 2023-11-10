# Create Capture History (CH) from CES data

require(reshape2)

# produce capture history
# NAs_in_missing_visits: 
#  - TRUE - put NAs instead of zero's where visits are missing. In this case, ces must contain ALL of the data
#    (all captures of all species) to figure out the visits (there was no visit without any bird caught).
#    You can always suply all data in ces, it doesn't speed up the process much to give a subset. Setting
#    this parameter to FALSE does speed it up significantly though.
#  - FALSE - just keep zero's whenever the birds wasn't caught (regardless whether there was a visit or not)
create_CH <- function (ces, filter, YRS, NAs_in_missing_visits = TRUE)
{

	ces$dummy <- TRUE

	# cast: drop = FALSE and factor(year,..) fills in all missing factor levels, if any
	CH0 <- acast(ces[filter & ces$year %in% YRS,], exID ~ factor(year, levels = YRS),
		any, value.var = "dummy", drop = FALSE)
	
	if (!NAs_in_missing_visits) { 
		CH <- ifelse(CH0, 1, 0)
		return(CH)
	}

# modify capture history to put NAs instead of zero's where visits are missing
	ces_yr_visits <- sqldf("
	select site_ces, year
	from ces
	group by site_ces, year
	")

	CH_visits0 <- sqldf("
	select exID, ces_yr_visits.year
	from CH0d join ces using (exID)
		join ces_yr_visits using (site_ces)
	group by exID, ces_yr_visits.year
	")
	CH_visits0$dummy <- TRUE
	CH_visits <- acast(CH_visits0[CH_visits0$year %in% YRS,], exID ~ factor(year, levels = YRS), any, value.var = "dummy", drop = FALSE)
	CH <- ifelse(CH0, 1, ifelse(CH_visits, 0, NA))

	CH
}


# Transience trick (sometimes called extended capture histories): 
# Behind the first occurrence, insert the known residency status (num.captures.first.year > 1)
extend_CH <- function (CH, num.captures.first.year)
{
	# first we will convert CH to strings, we will do it on the strings and then we convert it back
	ch_tmp0 <- apply(ifelse(is.na(CH), ".", as.character(CH)), 1, paste0, collapse = "")
	xx <- cbind(ch_tmp0, paste0("1", ifelse(num.captures.first.year > 1, "1", "0")))
	ch_tmp <- apply(xx, 1, function(x) sub("1", x[2], x[1]))
	stopifnot(all(nchar(ch_tmp) == ncol(CH)+1))
	require(stringr)
	CH2 <- str_split_fixed(ch_tmp, "", ncol(CH)+1)
	class(CH2) <- "numeric" # warning "NAs introduced by coercion" is allright here
	if (1) { # small test
		CH.test <- str_split_fixed(ch_tmp0, "", ncol(CH))
		stopifnot(all(CH == CH.test))
		rm(CH.test)
	}
	CH2 # output: CH2 - with the "transients trick" :)
}

# in how many visits was the individual captured during the occasion of the last capture (occasion here referred to as "year", but can be different time unit)
num_captures_last_year <- function (ces, filter, YRS)	
{
	CH0sum <- acast(ces[filter & ces$year %in% YRS,], exID ~ year,
		function (x) length(unique(x)), value.var = "visit")

	get.last <- function(x) max(which(x>0))
	apply(CH0sum, 1, function (x) x[get.last(x)])
}


# in how many visits was the individual captured during the occasion of the first capture (occasion here referred to as "year", but can be different time unit)
num_captures_first_year <- function (ces, filter, YRS)
{
	CH0sum <- acast(ces[filter & ces$year %in% YRS,], exID ~ year,
		function (x) length(unique(x)), value.var = "visit")

	get.first <- function(x) min(which(x>0))
	apply(CH0sum, 1, function (x) x[get.first(x)])
}









