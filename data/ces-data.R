# Load the data (individual capture records from CES - Constant Effort Sites)

require(sqldf)

load(file = "ces-data.Rdata")

ces$date <- as.Date(paste(ces$day, ces$month, ces$year, sep = "."), format="%d.%m.%Y")

ces$exID <- ces$ringID # exID: unique ID of an individual

ces_visits <- sqldf("
select site_ces, date, yday, year, count(*) as num_obs
from ces
group by site_ces, date")

