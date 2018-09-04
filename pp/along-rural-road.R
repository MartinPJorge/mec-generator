# @note: this script is used to generate the Hoces-del-Cabriel road antennas
# be sure you execute it next to this script location

source("gen-utils-clean.R")
library(ggmap)
library(spatstat)
library(graphics) 
library(rjson)
library(cluster)
library(latex2exp)
library(metR)

REGIONS <- "../data/regions.json"
REGION_NAME <- "Hoces-del-Cabriel"
REGION <- "../data/antennas/Hoces-del-Cabriel/Hoces-del-Cabriel.csv"
HOCES_DEL_CABRIEL_ROAD_DOWN_LAT <- 39.5
TELEFONICA_NET <- 7
RURAL_ANTENNAS_DIS <- 1000 # 1 antenna per square km


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
# Remove below road
regionAntennas <- subset(regionAntennas,
                         regionAntennas$lat > HOCES_DEL_CABRIEL_ROAD_DOWN_LAT)

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
                A3Markpoints = region_$A3Markpoints)
    break
  }
}



currLon <- region$bl$lon
# Order the region antennas longitude values
antennasOrderLons <- regionAntennas[order(regionAntennas$lon),]
orderedLons <- antennasOrderLons$lon
currLon <- destination(lat = antennasOrderLons[1,]$lat,
                                 lon = currLon, bearing = 90,
                                 distance = RURAL_ANTENNAS_DIS / 2)$lon2
positionedLons <- c(currLon)
while (currLon < region$br$lon) {
  currLon <- destination(lat = antennasOrderLons[1,]$lat,
                                   lon = currLon, bearing = 90,
                                   distance = RURAL_ANTENNAS_DIS)$lon2
  positionedLons <- c(positionedLons, currLon)
}

# Obtain the latitudes of the generated road antennas' longitudes
roadAntennas <- roadLatitudes(pointsLongitudes = positionedLons,
                                milestones = region$A3Markpoints)

# Append the generated road antennas to the existing ones
roadAntennasW <- owin(xrange = c(region$bl$lon - 1, region$tr$lon + 1),
                      yrange = c(region$bl$lat - 10, region$tr$lat + 10))
roadAntennasPPP <- ppp(x = roadAntennas$lon, y = roadAntennas$lat,
                       window = roadAntennasW)
roadAndAllAntennas <- appendSmallCells(antennasDf = allRegionAntennas,
                                 smallCellsPP = roadAntennasPPP,
                                 operatorNET = TELEFONICA_NET)

if (isScript) {
  write.csv(x = roadAndAllAntennas, file = OUT_CSV_PATH)
}


if (!isScript) {
  # Get the map region
  mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
                 right = region$tr$lon, top = region$tr$lat)
  map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
                maptype = region$plotDetails$mapType)
  
  # Get the mark locations for eye reference of the lines covering the road
  markLats <- c()
  markLons <- c()
  for (i in 1:length(region$A3Markpoints)) {
    markLats <- c(markLats, region$A3Markpoints[[i]]$lat)
    markLons <- c(markLons, region$A3Markpoints[[i]]$lon)
  }
  marks_ <- data.frame(lon = markLons, lat = markLats)
  
  # Plot the map
  ggmap(map) + 
    geom_point(data = roadAntennas, shape = 17, size = 1) +
    geom_point(data = marks_, shape = 13, size = 2) +
    geom_point(data = allRegionAntennas, shape = 18, size = 5)
}


