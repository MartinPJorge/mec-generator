source("gen-utils-clean.R")
library(rjson)

# Default global variable
ANTENNAS <- TRUE
PEOPLE <- FALSE
REGIONS <- "../data/regions.json"
REGION_CSV <- "../data/antennas/Madrid-center/Madrid-center.csv"
REGION_NAME <- "Madrid-center"
LON_N <- 100 # number of divisions in the longitude
LAT_N <- 100 # number of divisions in the latitude

# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 4) {
    stop(cat("Arguments to receive are: ",
      "antenna|people regionId longitudeN latitudeN"))
  } else if(args[1] != "antenna" & args[1] != "people") {
    stop("First argument must be 'antenna' or 'people'")
  } else {
    REGION_NAME <- args[2]
    LON_N <- as.numeric(args[3])
    LAT_N <- as.numeric(args[4])
    if(args[1] == "people") {
      ANTENNAS <- FALSE
      PEOPLE <- TRUE
      REGION_CSV <- list.files(path = paste("../data/people/", REGION_NAME,
                                            sep = ""))
    }
    REGION_CSV <- paste("../data/antennas/", REGION_NAME, "/",
                           REGION_NAME, ".csv", sep = "")
  }

# Load files
regions <- fromJSON(file = REGIONS)
regionAntennas <- read.csv(REGION_CSV)

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
                populations = region_$populations)
    break
  }
}


################ START HERE
outGriddedJson <- "../data/"
if (ANTENNAS) {
  outGriddedJson <- paste(outGriddedJson, "antennas/", 
    REGION_NAME, "/", 
    "lat", toString(LAT_N), "lon", toString(LON_N), ".json", sep = "")
  gridData(inCSV = REGION_CSV, antennas = TRUE, region = region,
               lonN = LON_N, latN = LAT_N, griddedJSON = outGriddedJson)
} else {
  for (i in 1:length(REGION_CSV)) {
    outGriddedJson <- paste(outGriddedJson, "people/",
      REGION_NAME,
      "/lat", toString(LAT_N), "lon", toString(LON_N), 
      "/sim", as.numeric(as.POSIXct(Sys.time())), "-", toString(i),
      ".json", sep = "")
    gridData(inCSV = REGION_CSV[[i]], antennas = TRUE, region = region,
                 lonN = LON_N, latN = LAT_N, griddedJSON = outGriddedJson)
  }
}


