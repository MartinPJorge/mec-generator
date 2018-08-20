source("gen-utils-clean.R")
library(rjson)

# Default global variable
ANTENNAS <- TRUE
PEOPLE <- FALSE
REGIONS <- "../data/regions.json"
DATA_CSV <- "../data/antennas/Madrid-center/Madrid-center.csv"
REGION_NAME <- "Madrid-center"
LON_N <- 100 # number of divisions in the longitude
LAT_N <- 100 # number of divisions in the latitude

# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 4) {
    stop(cat("Arguments to receive are: ",
      "antenna|people regionId longitudeN latitudeN\n"))
  } else if(args[1] != "antenna" & args[1] != "people") {
    stop("First argument must be 'antenna' or 'people\n'")
  } else {
    REGION_NAME <- args[2]
    LON_N <- as.numeric(args[3])
    LAT_N <- as.numeric(args[4])
    if(args[1] == "people") {
      ANTENNAS <- FALSE
      PEOPLE <- TRUE
      peopleDir <- paste("../data/people/", REGION_NAME, sep = "")
      DATA_CSV <- list.files(path = peopleDir, pattern = "*.csv",
                             full.names = TRUE)
      if (length(DATA_CSV) == 0) {
        stop(paste("No people CSVs under: ", peopleDir,
                   " please execute people-sandbox.R to generate them\n",
                   sep = ""))
      }
    }
    else
      DATA_CSV <- paste("../data/antennas/", REGION_NAME, "/",
                             REGION_NAME, ".csv", sep = "")
  }


# Find the selected region
regions <- fromJSON(file = REGIONS)
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
if (ANTENNAS) {
  outGriddedJson <- paste("../data/antennas/", 
    REGION_NAME, "/", 
    "griddedLon", toString(LON_N), "lat", toString(LAT_N), ".json", sep = "")
  print(paste("Gridding ", DATA_CSV, sep = ""))
  gridData(inCSV = DATA_CSV, antennas = TRUE, region = region,
               lonN = LON_N, latN = LAT_N, griddedJSON = outGriddedJson)
} else {
  for (i in 1:length(DATA_CSV)) {
    outGriddedJson <- paste("../data/people/",
      REGION_NAME,
      "/griddedLon", toString(LON_N), "lat", toString(LAT_N), 
      "sim", as.numeric(as.POSIXct(Sys.time())), "-", toString(i),
      ".json", sep = "")
    print(paste("Gridding ", DATA_CSV[[i]], sep = ""))
    gridData(inCSV = DATA_CSV[[i]], antennas = FALSE, region = region,
                 lonN = LON_N, latN = LAT_N, griddedJSON = outGriddedJson)
  }
}


