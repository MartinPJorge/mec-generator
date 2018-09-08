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
OUT_MAT <- NULL
OUT_LONS <- NULL
OUT_LATS <- NULL
CLI <- FALSE # flag to tell if file is executed from CLI



# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 6) {
    stop(paste("Arguments to receive are: ",
      "regionId longitudeSamp latitudeSamp matOUT lonsOUT latsOUT"))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    LON_SAMPLES <- as.numeric(args[2])
    LAT_SAMPLES <- as.numeric(args[3])
    OUT_MAT <- args[4]
    OUT_LONS <- args[5]
    OUT_LATS <- args[6]
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


# Generate the intensity matrix of people used for the lambda of people
peopleF <- genPopuManta(townCents = townCents, townRads = townRads,
                        townPopus = townPopus, townDisps = townDisps)
intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = LAT_SAMPLES,
                        longisLen = LON_SAMPLES, intF = peopleF) 

# Normalize the lambda matrix
#rePopMat <- log(intMat$matrix + 1)
rePopMat <- intMat$matrix
rePopMat <- rescaleRange(matrix = rePopMat, toA = 1, toB = 1.07)

# Aleatorize with small scale to prevent equal marks
minRan <- min(rePopMat) / 1000000
randMat <- replicate(nrow(rePopMat), runif(ncol(rePopMat), min = minRan,
                                           max = 2 * minRan))
rePopMat <- rePopMat + randMat


if (CLI) {
  write.csv(rePopMat, file = OUT_MAT)
  write(intMat$lonAxis, file = OUT_LONS)
  write(intMat$latAxis, file = OUT_LATS)
}

