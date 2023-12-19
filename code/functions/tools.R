# tools.R - various utility functions


# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
# from http://www.cookbook-r.com/Manipulating_data/Comparing_vectors_or_factors_with_NA/
# https://stackoverflow.com/a/61248968/684229
#
# !! pozor nehlida ruznou delku
compareNA <- function(v1,v2) {
	same <- (v1 == v2) | (is.na(v1) & is.na(v2))
	same[is.na(same)] <- FALSE
	return(same)
}

# same as compareNA, but with NULL for non-vector variables
compareNULL <- function(v1,v2) {
	stopifnot(length(v1) <= 1)
	stopifnot(length(v2) <= 1)
	(is.null(v1) && is.null(v2)) || (!is.null(v1) && !is.null(v2) && v1 == v2)
}

# mstart(), mstop(): pair of functions for measurement of elapsed time and amount of used memory
#
# bacha resetuje gc(), tj. vola gc(reset = TRUE)!
# absolute_time - obsolete way to print abs. time
# id != NULL voli nove rozhrani, ktere resi vnorene cally, matchuje mstop s mstart se stejnym id; tato varianta uz vraci i spravnou velikost 
#    pameti pri nested callech! 
#	id = NULL pro zpetnou kompatibilitu, kde se to bastli do jedne glob. promenne, a to pak nefunguje s vnorenymi cally ani pro pamet, ani pro cas!!
#	az zrusim moznost id = NULL pro zp. kompat., a hodim nove rozhrani i pro id = NULL, tak bude mozne mstart/mstop volat bez konkretniho id
#	a bude to fungovat v poho (id pro tu fcnost neni potreba, resi se to poradim na stacku; id se jen kontroluje zda sedi)
# mem_precise = TRUE: it will run actual garbage collection, which will make the maxMemMB value more precise (maxMemMBprocess should be good anyway).
#	the g.c. usually takes 0.2 sec, so it might not be suitable if a lot of mstart/mstop are called. PS: not sure if this is enough to make
#	maxMemMB precise, because the maxMemMBprocess value itself might be a bit higher if some memory was not gc-ed at the peak moment... just
#   food for thought, not sure about this.
#
# In the output, the tilde (~) will show if the value of maxMemMB is imprecise. 
# The difference between maxMemMBprocess - maxMemMB is exactly the amount of memory allocated before the matching mstart() call and this 
# is the thing which is subject to the mentioned (im)precision.
#
# !!! mem_precise = FALSE probably not working so well! See d:\tomas\acr\vypocet\results\gp-cv2\DOHRO_FINALE\models\PDP\pred-2D\pre2-chyba\logs\Carduelis_flammea_gp-outputs-log-LAPTOP-SAAK4JMO_2688.txt  :
# 		pred() took 2.2 sec, max memory used: 621.3Mb, whole process: 2576.5Mb.
# 		model_expand_predictions() took 0.1 sec, max memory used: 41.6Mb, whole process: 1996.5Mb.
# 		returning memory - gc() took 0.4 sec, max memory used: ~0.0Mb, whole process: 1955.1Mb.
# 		whole predict() took 3.3 sec, max memory used: 621.3Mb, whole process: 2576.5Mb.
# 		pred.PDP took 10119.3 sec, max memory used: ~351.6Mb, whole process: 3758.1Mb.
# - how can pred.PDP take only 351.6Mb of max memory used, when it is a calling function of the preceding ones?
#
# Careful when using these within tryCatch!! See test_mstart_exceptions*.R.
# Workaround for now: don't use mstart/mstop in tryCatch; or in case of redis workers etc., only use mstart/mstop within tryCatch.
#
# !!! memory reporting will not work if anyone will be calling gc(reset = TRUE) after mstart() (mozna udelat wrapper pro gc(), kterej bude pak
# updatovat max_mem na stacku? To by bylo robustni reseni!)
#
# TIP: if you want really good debugging, call this:
# mstartOptions(mem_precise = TRUE, report_id = TRUE)
# this will change these options for all calls of mstart/mstop anywhere!! :)
mstart <- function(absolute_time = FALSE, id = NULL, bw.compat = FALSE, quiet = FALSE, mem_precise = mstartOptions("mem_precise"))
{
	gc_max <- gc(full = FALSE) # get the maximum before reset
	if (mem_precise)
		gc1 <- gc(reset = TRUE) # full = TRUE => g.c. action will take place and eat around 0.2sec...  at least here it's necessary to free the memory, if we want maxMemMB to be real. maxMemMBprocess would be good anyway even if we didn't free the g.c. memory here.
	else
		gc1 <- gc(reset = TRUE, full = FALSE)
	start_time <- Sys.time()
	if (absolute_time && !quiet)
		print(start_time)
	t1 <- proc.time()
	if (bw.compat) { # is.null(id)) {
		# nejsem si jist, jestli jeste nejaky kod nespoleha na to, ze se to dava do glob. promennych. Nechavama to tu pro zp. kompat.
		gc1 <<- gc1
		t1 <<- t1
		if (absolute_time)
			start_time <<- start_time
	} else {
		# nove rozhrani bez glob promennych. Ani by nepotrebovalo to id pro identifikaci: staci pouzivat stack a pririzenou "nestedness". id je spis pro kontrolu
		# uz jen a proto, abych specifikoval ze chci nove rozhrani. V dalsi verzi se muze toto udelat jako default a na zp. kompat. se vykaslat.
		if (!exists("mstart_mstop_data"))
			mstart_mstop_data <<- mstartDataInit()
		mstart_mstop_data <<- within(mstart_mstop_data, {
			if (i > 0) {
				#cat("i = ", i, ": going to do stack[[i]]$max_mem <- max(stack[[i]]$max_mem == ", stack[[i]]$max_mem, ", sum(gc_max[,6]) == ", sum(gc_max[,6]), ")\n")
				stack[[i]]$max_mem <- max(stack[[i]]$max_mem, sum(gc_max[,6])) # maximum memory eaten by the whole process since the opening mstart() (which this record corresponds to) up to the last mstart()
			}
			i <- i + 1
			stack[[i]] <- list(
				max_mem = 0,
				gc1 = gc1,
				t1 = t1,
				id = id,
				start_time = start_time,
				mem_precise = mem_precise
			)
		})
	}
} # tested, viz test_mstart_mstop-nested_calls,memory.R

mstartDataInit <- function () list(i = 0, stack = list(), options = list(mem_precise = FALSE, report_id = FALSE))

mstartOptions <- function (getOption = NULL, mem_precise = NULL, report_id = NULL)
{
	if (!exists("mstart_mstop_data"))
		mstart_mstop_data <<- mstartDataInit()
	if (!is.null(getOption))
		return(mstart_mstop_data$options[[getOption]])
	if (!is.null(mem_precise))
		mstart_mstop_data$options$mem_precise <<- mem_precise
	if (!is.null(report_id))
		mstart_mstop_data$options$report_id <<- report_id		
}

mstop <- function(absolute_time = FALSE, id = NULL, bw.compat = FALSE, quiet = FALSE, report_id = mstartOptions("report_id")) 
{
	t2 <- proc.time()
	gc2 <- gc(full = FALSE)
	max_mem <- 0
	end_time <- Sys.time()
	mem_precise <- NA
	if (absolute_time && !quiet)
		print(end_time)
	if (!bw.compat) { # !is.null(id)) { # see comment in mstart
		if (!exists("mstart_mstop_data"))
			stop("Matching mstart() call not found")
		t1 <- NA
		gc1 <- NA
		start_time <- NA
		mstart_mstop_data <<- within(mstart_mstop_data, {
			if (length(stack) < 1)
				stop("mstop(): there is nothing on the stack - did you call mstart()?")
			if (!compareNULL(stack[[i]]$id, id))
				stop("mstop(): id = ", id, " doesn't match the one on the stack (" ,stack[[i]]$id, ") - non-matched mstart/mstop calls")
			t1 <<- stack[[i]]$t1
			gc1 <<- stack[[i]]$gc1
			max_mem <<- stack[[i]]$max_mem
			start_time <<- stack[[i]]$start_time
			mem_precise <<- stack[[i]]$mem_precise
			stack[[i]] <- NULL
			i <- i - 1
			if (i > 0) 
				stack[[i]]$max_mem <- max(stack[[i]]$max_mem, max_mem) # predej maximum o level vys!
			
		})
	}
	#print(t2 - t1)
	timeSec <- (t2 - t1)[3]
	maxMemMBprocess <- max(sum(gc2[,6]), max_mem)  # maximum memory eaten by the whole process since the matching mstart() call
		# max_mem = maximum memory eaten by the whole process since the matching mstart() up to the previous mstart() call
		# sum(gc2[,6]) = maximum memory eaten by the whole process since last mstart() call
	maxMemMB <- maxMemMBprocess - sum(gc1[,2]) # maximum memory used just by the code running since the matching mstart() call
		# sum(gc1[,2]) = memory used by the process before the matching mstart() call. If mem_precise = FALSE, then this be higher since it might 
		# include unused memory which hasn't been freed, because when mem_precise = FALSE, garbage collection doesn't take place in mstart().
	if (!quiet) {
		if (report_id && !is.null(id))
			cat(sprintf("id=\"%s\": ", id))
		cat(sprintf("%.1f sec, ", timeSec))
		if (mem_precise)
			cat(sprintf("max memory used: %.1fMb", maxMemMB))
		else
			cat(sprintf("max memory used: ~%.1fMb", maxMemMB))
		cat(sprintf(", whole process: %.1fMb.\n", maxMemMBprocess))
	}
	if (maxMemMB < 0) {
		cat("mstop debug: sum(gc2[,6]) = ", sum(gc2[,6]), ", max_mem = ", max_mem, ", max_mem_process = ", maxMemMBprocess, "\n")
		cat("gc2:\n")
		print(gc2)
		cat("gc1:\n")
		print(gc1)
		if (0) {
			hostname <- Sys.info()["nodename"]
			dump.fn <- paste0("mstop_debug,", hostname, "_", Sys.getpid(), "_dump")
			warning("mstop debug: maxMemMB < 0: dumping stack, variables etc. to ", dump.fn, "\n")
			dump.frames(dumpto = dump.fn, to.file = TRUE, include.GlobalEnv = TRUE) # use debugger(loadVar("*", dump.fn)) to debug...
		}
	}	
	if (bw.compat) { # (is.null(id)) { # pro zp. kompat. radeji (projit vsechen kod co pouziva mstart/mstop a podivat se, jestli to nejde zrusit)
		gc2 <<- gc2
		t2 <<- t2
		if (absolute_time)
			end_time <<- end_time
	} else {
		#cat("I am here")
		return(invisible(list(timeSec = timeSec, maxMemMB = maxMemMB, maxMemMB_precise = mem_precise, maxMemMBprocess = maxMemMBprocess, startTime = start_time, endTime = end_time)))
	}
}


