source("gen-utils-clean.R")

# Default global variable
ANTENNAS <- TRUE
PEOPLE <- FALSE
REGIONS <- "../data/regions.json"
REGION_ANTENNAS <- "../data/antennas/Madrid-center/Madrid-center.csv"
REGION_NAME <- "Madrid-center"
LON_N <- 100 # number of divisions in the longitude
LAT_N <- 100 # number of divisions in the latitude

# Parse arguments if existing to change default global variables
# if(length(args) > 0)
#   if(length(args) != 4) {
#     stop(cat("Arguments to receive are: ",
#       "antenna|people regionId longitudeN latitudeN"))
#   } else if(args[1] != "antenna" & args[1] != "people") {
#     stop("First argument must be 'antenna' or 'people'")
#   } else {
#     if(args[1] == "people") {
#       ANTENNAS <- FALSE
#       PEOPLE <- TRUE
#     }
#     REGION_NAME <- args[2]
#     REGION_ANTENNAS <- paste("../data/antennas/", REGION_NAME, "/",
#                            REGION_NAME, ".csv", sep = "")
#     LON_N <- args[3]
#     LAT_N <- args[4]
#   }


# Load files
regions <- fromJSON(file = REGIONS)
regionAntennas <- read.csv(REGION_ANTENNAS)

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
  gridAntennas(antennasCSV = REGION_ANTENNAS, region = region,
               lonN = LON_N, latN = LAT_N, griddedJSON = outGriddedJson)
} else {
  outGriddedJson <- paste(outGriddedJson + "people/" +
    REGION_NAME + 
    "/lat", toString(LAT_N), "lon", toString(LON_N), 
    "/sim", as.numeric(as.POSIXct(Sys.time())), ".json", sep = "")
}


