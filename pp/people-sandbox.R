source("gen-utils-clean.R")
library(pracma)
library(SDMTools)
library(ggmap)
library(spatstat)
library(graphics) 
library(rjson)
library(cluster)
library(latex2exp)
library(metR)

REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
N_SIMS <- 5 # number of simulations to generate people in selected region
LON_SAMPLES <- 100
LAT_SAMPLES <- 100
CLI <- FALSE # flag to tell if file is executed from CLI


# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 4) {
    stop(paste("Arguments to receive are: ",
      "regionId numSims longitudeSamp latitudeSamp"))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    N_SIMS <- as.numeric(args[2])
    LON_SAMPLES <- as.numeric(args[3])
    LAT_SAMPLES <- as.numeric(args[4])
  }


# Load files
regions <- fromJSON(file = REGIONS)

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

# Obtain the region population areas and its details
townPopus <- c()
townRads <- c()
townDisps <- c()
popuLons <- c()
popuLats <- c()
for (population in region$populations) {
  townPopus <- c(townPopus, population$population)
  townRads <- c(townRads, population$radius)
  coords <- getCoords(population$center)
  popuLons <- c(popuLons, coords$lon)
  popuLats <- c(popuLats, coords$lat)
  townDisps <- c(townDisps, population$dispersion)
}
townCents <- data.frame(lat = popuLats, lon = popuLons)


##################### DATA GENERATION #####################
# Generate the intensity function
peopleF <- genPopuManta(townCents = townCents, townRads = townRads,
                        townPopus = townPopus, townDisps = townDisps)
intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = LAT_SAMPLES,
                        longisLen = LON_SAMPLES, intF = peopleF) 
intFrame <- intMat2Frame(intMat$matrix, intMat$latAxis, intMat$lonAxis)
peopleIntpF <- genInterpInt(intMat$matrix, intMat$lonAxis, intMat$latAxis)

# Generate people
ppWindow <- owin(xrange = c(region$bl$lon, region$tr$lon),
                   yrange = c(region$bl$lat, region$tr$lat))
print(paste("Generating", toString(N_SIMS), "simulations of points..."))
peoplePoints <- rpoispp(lambda = peopleIntpF, win = ppWindow, nsim = N_SIMS)
for (i in 1:N_SIMS) {
  peopleOutF <- paste("../data/people/", REGION_NAME,
                      "/lon", toString(LON_SAMPLES),
                      "lat", toString(LAT_SAMPLES), "-",
                      toString(as.numeric(as.POSIXct(Sys.time()))),
                      "-", toString(i), ".csv", sep = "")
  print(paste("Simulation ", toString(i), " stored at ", peopleOutF, sep = ""))
  popuToCSV(popuPoints = peoplePoints[[i]], outCSV = peopleOutF)
}
##########################################################


# PLOT if not run from CLI
if (!CLI) {
  # Get the logarithm of intensities so small villeages are shown
  logIntens <- c()
  for (row in 1:nrow(intFrame)) {
    v <- intFrame$intensity[row]
    if (v != 0) {
      v <- log(v)
    }
    logIntens <- c(logIntens, v)
  }
  loggedIntF <- data.frame(lon = intFrame$lon, lat = intFrame$lat,
                           intensity = logIntens)
  
  mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
                 right = region$tr$lon, top = region$tr$lat)
  map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
                maptype = region$plotDetails$mapType)
  
  ggmap(map) + 
    geom_contour_fill(data = loggedIntF, aes(z=intensity), alpha = 0.15) +
    scale_fill_gradient(name = TeX("$\\lambda (u)$"), low = "gray", high = "black") +
    scale_colour_gradient(low = "gray", high = "black") +
    geom_point(data = data.frame(lon = people$x, lat = people$y))
}



        