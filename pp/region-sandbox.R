source("gen-utils-clean.R")
library(ggmap)
library(spatstat)
library(graphics) 
library(rjson)
library(cluster)
library(latex2exp)

REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
REGION <- "../data/Madrid-center.csv"
# REGION_NAME <- "Cobo-Calleja"
# REGION <- "../data/Cobo-Calleja.csv"

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
                tr = getCoords(region_$tr))
    break
  }
}

# Extract the location of the antennasd
HEAD_ANTENNAS <- 10 # <- vary this to know with how many antennas can work
antennasLoc <- head(regionAntennas, HEAD_ANTENNAS)

# Gener
intensityF <- gen_manta(antennasLoc, dgamma, shape=100, scale=10)
intFrame <- genIntFrame(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = 50, longisLen = 50, intF = intensityF) 



# Draw the map and the 
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = 13, source = "stamen",
               maptype = "toner-lite")
ggmap(map) + 
  stat_contour(data=intFrame, aes(z=intensity, fill=intensity, color=..level..),
               geom='polygon', alpha=0.2) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  labs(color = TeX("$\\lambda (u)$")) +
  geom_point(data = antennasLoc)
#  stat_contour(data=intFrame, aes(z=intensity)) +
#  scale_fill_gradient(low = "gray", high = "black")



ggmap(map)


qmplot(lon, lat, data = antennasLoc, maptype = "toner-lite")
