source("gen-utils-clean.R")
library(rjson)

REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
REGION <- "../data/antennas/Madrid-center/Madrid-center.csv"
TELEFONICA_NET <- 7


# Parse arguments if existing to change default global variables
isScript <- FALSE
OUT_CSV_PATH <- NULL
args <- commandArgs(trailingOnly=TRUE)
print(length(args))
if(length(args) > 0)
  if(length(args) != 1) {
    stop(cat("Arguments to receive are: ",
      "OUT_CSV_PATH <- the CSV whith the generated antennas appended\n"))
  } else {
    OUT_CSV_PATH <- args[1]
    isScript <- TRUE
  }


# Load files
regions <- fromJSON(file = REGIONS)
regionAntennas <- read.csv(REGION)
regionAntennas <- subset(regionAntennas, regionAntennas$radio == "LTE")
allRegionAntennas <- subset(regionAntennas, regionAntennas$net == TELEFONICA_NET)


# Find the selected region
region <- NULL
for (region_ in regions$regions) {
  if (region_$id == REGION_NAME) {
    region <- list(id = region_$id,
                bl = getCoords(region_$bl),
                br = getCoords(region_$br),
                tl = getCoords(region_$tl),
                tr = getCoords(region_$tr),
                repulsionRadius = region_$repulsionRadius,
                plotDetails = region_$plotDetails,
                skycrapers = region_$skycrapers)
    break
  }
}


# Append each femto-cell to the antennas CSV
femtoLons <- c()
femtoLats <- c()
for (s in 1:length(region$skycrapers)) {
  skycraper <- region$skycrapers[[s]]
  femtoCells <- 4 * skycraper$floors
  femtoLons <- c(femtoLons, rep(skycraper$lon, femtoCells))
  femtoLats <- c(femtoLats, rep(skycraper$lat, femtoCells))
}
antennasAndFemtos <- manualAntennaAppend(antennas = regionAntennas,
                                         newAntennas = data.frame(
                                           lat = femtoLats, lon = femtoLons),
                                         newRadio = "femto-cell", 
                                         newNetOperator = TELEFONICA_NET)

if (isScript) {
  write.csv(antennasAndFemtos, file = OUT_CSV_PATH)
} else {
  View(antennasAndFemtos)
}
