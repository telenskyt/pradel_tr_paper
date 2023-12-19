
library(reshape2)

# @description 
# Create capture history (CH) from individual capture data.
#
# @param cap - table (data.frame) with individual captures (might be e.g. all captures of certain species) - one line per capture. 
#		Multiple captures allowed. Following columns are assumed:
#		- exID - unique ID of the individual animal
#		- time - integer, the temporal occasion (e.g. a year); can start at any number; +1 means next occasion.
#
# @param occasions a vector of temporal occasions (values of the `time` column) for which to present the capture histories. This is useful to specify
#		explicitly e.g. in case the `cap` is just a selection (e.g. for some site), but you want a wider range of occasions to be presented than covered by this site.
# 
# @param visits - used when `NAs_in_missing_visits` is TRUE (not implemented yet). It is a data.frame with the columns `site` and `time` (same as in the `cap` parameter), 
#		one row for each site & time combination where the visit (capture occasion) was performed.
#
# @param NAs_in_missing_visits
#  - TRUE - (not implemented yet) put NAs instead of zero's where visits are missing. In this case, parameter `visits` must be supplied.
#  - FALSE - just keep zero's whenever the animal wasn't caught (regardless whether there was a visit or not)
create_CH <- function (cap, occasions = min(cap$time):max(cap$time), visits = NULL, NAs_in_missing_visits = TRUE)
{
	cap$dummy <- TRUE

	# cast: drop = FALSE and factor(time,..) fills in all missing factor levels, if any
	CH0 <- acast(cap, exID ~ factor(time, levels = occasions), any, value.var = "dummy", drop = FALSE)
	
	if (!NAs_in_missing_visits) { 
		CH <- ifelse(CH0, 1, 0)
		return(CH)
	}
	
	stop("The variant of NAs_in_missing_visits = TRUE is not implemented yet.")
		# Internal comment: it is implemented in my version of d:\tomas\ces\ces-data-capture_history.R, however this code was not adapted to this new interface yet.
	stopifnot(!is.null(visits))
}

# @description
# Transience trick (sometimes called extended capture histories): 
# For each individual, behind the first capture (first `1` in capture history), insert the residency status known during the occasion of first capture (`1` if known.resid.first.cap is TRUE, `0` otherwise) 
#
# @param CH capture history - a matrix - rows ~ individuals, columns ~ temporal occasions
#
# @param known.resid.first.cap a boolean vector that for each individual denotes whether 
# 	it was a confirmed resident during the occasion of the first capture (e.g. by capturing it during multiple sub-occasions within the occasion of the first capture, or by direct cues confirming residency)
#		(for more details see `Model extension 2: Confirmation of residency status within the first capture occasion` in the paper
#		 Telenský, T., Storch, D., Klvaňa, P., Reif, J. Extension of Pradel capture-recapture survival-recruitment model accounting for transients. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14262)
#
# @details Rows in both arguments - CH and num.captures.first.year - must correspond!
#
# @return Returns the extended capture history
extend_CH <- function (CH, known.resid.first.cap)
{
	stopifnot(!is.null(known.resid.first.cap) && all(!is.na(known.resid.first.cap)))
		# the function could be easily adapted to allow NAs also in this parameter, as it allows in CH, but it is not done at the moment
		
	# first we will convert CH to strings, we will do it on the strings and then we convert it back
	ch_tmp0 <- CH_to_char(CH)
	xx <- cbind(ch_tmp0, paste0("1", ifelse(known.resid.first.cap, "1", "0")))
	ch_tmp <- apply(xx, 1, function(x) sub("1", x[2], x[1]))
	stopifnot(all(nchar(ch_tmp) == ncol(CH)+1))	
	CH2 <- char_to_CH(ch_tmp)

	if (1) { # small test
		CH.test <- char_to_CH(ch_tmp0)
		stopifnot(all(CH == CH.test))
		rm(CH.test)
	}
	CH2 # output: CH2 - with the "transients trick" :)
}

# @description
# Convert capture histories (CH) to vector of character strings
#
# @param CH capture history - a matrix - rows ~ individuals, columns ~ temporal occasions
#
# @return vector of character strings; one character string for each individual (row of the matrix `CH`), 
# containing characters `0` (not captured), `1` (captured), and `.` (for NAs in capture history; typically means
# there was no visit - no monitoring, no capturing event).

CH_to_char <- function (CH)
{
	apply(ifelse(is.na(CH), ".", as.character(CH)), 1, paste0, collapse = "")
}

# @description
# Convert capture histories in the form of character strings to matrix form
#
# @param char_CH a vector of character strings, each character string is a capture history, 
# containing characters `0` (not captured), `1` (captured); !!! NOTE THAT IT DOESN'T SUPPORT `.` NOW YET!!(`.` typically means there was no visit - no monitoring, no capturing event).
#
# @return capture history in the form of matrix, with 0s and 1s (NAs not supported yet)
char_to_CH <- function (char_CH)
{
	library(stringr)
	stopifnot(length(char_CH) > 0)
	# test the assumption that all the string capture histories are of the same length:
	n <- nchar(char_CH[1])
	stopifnot(all(nchar(char_CH) == n))
	CH2 <- str_split_fixed(char_CH, "", n)
	class(CH2) <- "numeric" # warning "NAs introduced by coercion" is allright here
	CH2
}
