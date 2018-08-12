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
                tr = getCoords(region_$tr),
                repulsionRadius = region_$repulsionRadius)
    break
  }
}

# Extract the location of the antennasd
HEAD_ANTENNAS <- 100 # <- vary this to know with how many antennas can work
antennasLoc <- head(regionAntennas, HEAD_ANTENNAS)

# Intensity function generation
# intensityF <- gen_manta(antennasLoc, factor=1e20, dgamma, shape=1, scale=10000)
intensityF <- gen_manta(centers = antennasLoc, factor = 1000,
                        centers_intensity = stepIntensity,
                        a = 100, b = 150, c = 650, d = 700)
timestamp()
intFrame <- genIntFrame(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = 50, longisLen = 50, intF = intensityF) 
timestamp()

# Poisson point process generation for the PoPs
ppWindow <- owin(xrange = c(region$bl$lon, region$tr$lon),
                   yrange = c(region$bl$lat, region$tr$lat))
timestamp()
pops <- rpoispp(intensityF, win = ppWindow)
timestamp()
rThin <- region$repulsionRadius
timestamp()
pops_th <- rthin(X = pops, P = rad_thin, r = rThin, nsim = 1, drop=TRUE)
 

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
  geom_point(data = antennasLoc, aes(shape="triangle"))
#  geom_point(data = data.frame(lat = pops_th$y, lon = pops_th$x),
#             aes(shape="circle"))

#  stat_contour(data=intFrame, aes(z=intensity)) +
#  scale_fill_gradient(low = "gray", high = "black")


######### TEMPTATIVE FIX https://github.com/eliocamp/metR/blob/master/R/stat_contour_fill.R#L36
ggmap(map) + 
  stat_contour_fill(data=intFrame, aes(z=intensity, fill=intensity, color=..level..),
               geom='polygon', alpha=0.2) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  labs(color = TeX("$\\lambda (u)$")) +
  geom_point(data = antennasLoc, aes(shape="triangle"))




########## NOW WITH OUR DATA based in https://github.com/hypertidy/contourPolys
intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                        longis = c(region$bl$lon, region$tr$lon),
                        latisLen = 50, longisLen = 50, intF = intensityF) 
z <- intMat$matrix
y <- intMat$latAxis
x <- intMat$lonAxis

levels <- pretty(range(z), n = 20)
p <- contourPolys::fcontour(x, y, z, levels)
m <- cbind(lon = unlist(p[[1]]), 
           lat = unlist(p[[2]]), 
           lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
           upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
           g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 

gd <- as.data.frame(m)

ggmap(map) +
  geom_polygon(data = gd, aes(x = lon, y = lat,
                              group = g, fill = upper, alpha = 1)) +
  scale_fill_gradient(low = "gray", high = "black") +
  geom_point(data = antennasLoc, aes(shape="b")) +
  labs(fill = TeX("$\\lambda (u)$"))
  

  # ggmap(map) + 
#   stat_contour(data=intFrame, aes(z=intensity, fill=intensity, color=..level..),
#                geom='polygon', alpha=0.2) +
#   scale_fill_gradient(low = "gray", high = "black") +
#   scale_colour_gradient(low = "gray", high = "black") +
#   labs(color = TeX("$\\lambda (u)$")) +
#   geom_point(data = antennasLoc, aes(shape="triangle")) +
#   geom_point(data = data.frame(lat = pops_th$y, lon = pops_th$x),
#              aes(shape="circle"))
