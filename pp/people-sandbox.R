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
REGION <- "../data/Madrid-center/Madrid-center.csv"
# REGION_NAME <- "Cobo-Calleja"
# REGION <- "../data/Cobo-Calleja/Cobo-Calleja.csv"
# REGION_NAME <- "Valle-Ayora-Cofrentes"
# REGION <- "../data/Valle-Ayora-Cofrentes/Valle-Ayora-Cofrentes.csv"

# Load files
regions <- fromJSON(file = REGIONS)
regionAntennas <- read.csv(REGION)

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


# Generate the intensity function
peopleF <- genPopuManta(townCents = townCents, townRads = townRads,
                        townPopus = townPopus, townDisps = townDisps)

latitudeSamples <- 100
longitudeSamples <- 100
intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = latitudeSamples,
                        longisLen = longitudeSamples, intF = peopleF) 
intFrame <- intMat2Frame(intMat$matrix, intMat$latAxis, intMat$lonAxis)
peopleIntpF <- genInterpInt(intMat$matrix, intMat$lonAxis, intMat$latAxis)
avgPeopleValle <- inhIntM(lon1 = region$bl$lon, lat1 = region$bl$lat,
                     lon2 = region$tr$lon, lat2 = region$tr$lat,
                     intensityF = peopleF)

# Generate people
ppWindow <- owin(xrange = c(region$bl$lon, region$tr$lon),
                   yrange = c(region$bl$lat, region$tr$lat))
people <- rpoispp(lambda = peopleIntpF, win = ppWindow, nsim = 1)


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



        