source("gen-utils-clean.R")
library(ggmap)
library(spatstat)
library(graphics) 
library(rjson)
library(cluster)
library(latex2exp)
library(metR)

REGIONS <- "../data/regions.json"
# REGION_NAME <- "Madrid-center"
# REGION <- "../data/antennas/Madrid-center/Madrid-center.csv"
# REGION_NAME <- "Cobo-Calleja"
# REGION <- "../data/antennas/Cobo-Calleja/Cobo-Calleja.csv"
# REGION_NAME <- "Valle-Ayora-Cofrentes"
# REGION <- "../data/antennas/Valle-Ayora-Cofrentes/Valle-Ayora-Cofrentes.csv"
REGION_NAME <- "Hoces-del-Cabriel"
REGION <- "../data/antennas/Hoces-del-Cabriel/Hoces-del-Cabriel.csv"

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
                plotDetails = region_$plotDetails)
    break
  }
}

# Extract the location of the antennasd
HEAD_ANTENNAS <- 100 #dim(regionAntennas)[1] # <- vary this to know with how many antennas can work
antennasLoc <- head(regionAntennas, HEAD_ANTENNAS)
# if (REGION_NAME == "Madrid-center") {
#   antennasLoc <- subset(antennasLoc, antennasLoc$radio == "LTE")
# }

# Intensity function generation
# intensityF <- gen_manta(antennasLoc, factor=1e20, dgamma, shape=1, scale=10000)
repRad <- region$repulsionRadius
intensityF <- gen_manta(centers = antennasLoc, factor = 1e5,
                        centers_intensity = stepIntensity,
                        a = repRad, b = repRad + 50,
                        c = repRad + 550, d = repRad + 600)

# Generate points with custom MaternII
latitudeSamples <- 50
longitudeSamples <- 50
timestamp()
intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = latitudeSamples,
                        longisLen = longitudeSamples, intF = intensityF) 
timestamp()
intFrame <- intMat2Frame(intMat$matrix, intMat$latAxis, intMat$lonAxis)
lambIntp <- genInterpInt(intMat = intMat$matrix, lonAxis = intMat$lonAxis,
             latAxis = intMat$latAxis)
val <- lambIntp(-3.689918, 40.42288)
rThin <- 1000
ppWindow <- owin(xrange = c(region$bl$lon, region$tr$lon),
                   yrange = c(region$bl$lat, region$tr$lat))
timestamp()
pops_th <- jorgeMaternIImapI(lambda = lambIntp, win = ppWindow,  r = rThin)
timestamp()

# for (i in 1:pops_th$n) {
#   sprintf("(%f, %f)", pops_th[i]$x, pops_th[i]$y)
#   for (j in 1:pops_th$n) {
#     if(j != i) {
#       print(distanceMulti(lat1 = pops_th$y[i], lon1 = pops_th$x[i],
#                     lat2 = pops_th$y[j], lon2 = pops_th$x[j]))
#     }
#   }
#   print("")
# }

# ##### HARD-CORE MATERNI VERSION #####
# # Generate the samples to draw the intensity
# intFrame <- genIntFrame(latis = c(region$bl$lat, region$tr$lat),
#                         longis = c(region$bl$lon, region$tr$lon),
#                         latisLen = latitudeSamples,
#                         longisLen = longitudeSamples, intF = intensityF) 
# # Poisson point process generation for the PoPs
# ppWindow <- owin(xrange = c(region$bl$lon, region$tr$lon),
#                    yrange = c(region$bl$lat, region$tr$lat))
# pops <- rpoispp(intensityF, win = ppWindow)
# rThin <- region$repulsionRadius
# pops_th <- rthin(X = pops, P = rad_thin, r = rThin, nsim = 1, drop=TRUE)
#  
# 
# # Draw the map and the 
# mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
#                right = region$tr$lon, top = region$tr$lat)
# map <- get_map(mapRegion, zoom = 13, source = "stamen",
#               maptype = "toner-lite")
# #####################################

######### TEMPTATIVE FIX https://github.com/eliocamp/metR/blob/master/R/stat_contour_fill.R#L36
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
              maptype = region$plotDetails$mapType)

ggmap(map) + 
  stat_contour_fill(data=intFrame, aes(z=intensity, fill=intensity, color=..level..),
               geom='polygon', alpha=0.2) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  labs(color = TeX("$\\lambda (u)$")) +
  geom_point(data = antennasLoc, aes(shape="triangle")) +
  geom_point(data = data.frame(lat = pops_th$y, lon = pops_th$x),
             aes(shape="circle"))


map <- get_map(mapRegion, zoom = 13, source = "stamen",
              maptype = "toner-lite")
ggmap(map) + 
  geom_point(data = antennasLoc, shape = 17, size = 3)


