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
MEC_GEN_PARAMS <- "../data/mec-gen-params.json"
LAT_SAMPLES <- 100
LON_SAMPLES <- 100
NUM_MECs <- 10
MACRO_CELLS_CSV <- "../data/antennas/Madrid-center/macro-cells.csv"
ANTENNAS_CSV <- "../data/antennas/Madrid-center/TEL-and-femtos.csv" # femto cells and LTE
OPERATOR <- 7 # net id of Telefonica
MEC_LOCATIONS_CSV <- "../data/antennas/Madrid-center/MEC-locations.csv"
CLI <- FALSE # flag to tell if file is executed from CLI



# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 9) {
    stop(paste("Arguments to receive are: ",
      "regionID latSamples lonSamples mecParamsJSON numMECs macro-cells.csv LTEandFEMTO.csv operatorNET MEClocations.csv"))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    MEC_GEN_PARAMS <- args[2]
    LAT_SAMPLES <- args[3]
    LON_SAMPLES <- arg[4]
    NUM_MECs <- as.numeric(args[5])
    MACRO_CELLS_CSV <- args[6]
    ANTENNAS_CSV <- args[7]
    OPERATOR <- as.numeric(args[8])
    MEC_LOCATIONS_CSV <- args[9]
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

# Read MEC generation parameters
mecParams <- fromJSON(file = MEC_GEN_PARAMS)

# Generate a data.frame with all the antennas
antennas <- read.csv(file = ANTENNAS_CSV)
antennas <- subset(antennas, antennas$net == OPERATOR)
macroCells <- read.csv(file = MACRO_CELLS_CSV)
lons <- c()
lats <- c()
radio <- c()
net <- c()
# Loop LTE and femto-cell antennas
for (row in 1:nrow(antennas)) {
  antenna <- antennas[row,]
  lons <- c(lons, antenna$lon)
  lats <- c(lats, antenna$lat)
  radio <- c(radio, as.character(antenna$radio))
  net <- c(net, antenna$net)
}
for (row in 1:nrow(macroCells)) {
  macroCell <- macroCells[row,]
  lons <- c(lons, macroCell$lon)
  lats <- c(lats, macroCell$lat)
  radio <- c(radio, as.character(macroCell$radio))
  net <- c(net, OPERATOR)
}
allAntennas <- data.frame(radio = radio, net = net, lon = lons, lat = lats)

# Generate the MEC server locations intensity function
maxDiss <- list(LTE = mecParams$LTE$maxDistance,
                macroCell = mecParams$smallCell$maxDistance,
                femtoCell = mecParams$femtoCell$maxDistance)

mecIntMat <- mecIntManhattan(lonL = region$bl$lon, lonR = region$tr$lon,
                latB = region$bl$lat, latT = region$tr$lat,
                lonSamples = LON_SAMPLES, latSamples = LAT_SAMPLES,
                antLons = allAntennas$lon, antLats = allAntennas$lat,
                antRadios = as.character(allAntennas$radio), maxDiss = maxDiss)
loggedMecInt <- intMat2Frame(intMat = mecIntMat$matrix,
                         latAxis = mecIntMat$latAxis,
                         lonAxis = mecIntMat$lonAxis)

# Get the map
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
              maptype = region$plotDetails$mapType)

ggmap(map) + 
  ggplot2::stat_contour(data = loggedMecInt, aes(z=intensity, fill=..level..),
               geom = "polygon", alpha=0.3) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  labs(fill = TeX("$log(1 + \\lambda (lon, lat))$")) 


# Obtain the MEC server locations using the antenna-decrease approach
numMECs <- 30
mecLocs <- mecLocationAntenna(mecIntMatrix = mecIntMat$matrix,
                              lonAxis = mecIntMat$lonAxis,
                              latAxis = mecIntMat$latAxis, maxDiss = maxDiss,
                              numMECs = NULL, antLons = allAntennas$lon,
                              antLats = allAntennas$lat,
                              antRadios = as.character(allAntennas$radio),
                              lonL = region$bl$lon, lonR = region$tr$lon,
                              latB = region$bl$lat, latT = region$tr$lat) 
cdfAntennaRed <- ecdf(mecLocs$pos$coveredAs)

mecLocs <- mecLocationDig(mecIntMatrix = mecIntMat$matrix,
                              lonAxis = mecIntMat$lonAxis,
                              latAxis = mecIntMat$latAxis, maxDiss = maxDiss,
                              numMECs = NULL, antLons = allAntennas$lon,
                              antLats = allAntennas$lat,
                              antRadios = as.character(allAntennas$radio),
                              lonL = region$bl$lon, lonR = region$tr$lon,
                              latB = region$bl$lat, latT = region$tr$lat) 
cdfDig <- ecdf(mecLocs$pos$coveredAs)

loggedMecLocsMat <- intMat2Frame(intMat = log(1 + mecLocs$modMat),
                         latAxis = mecIntMat$latAxis,
                         lonAxis = mecIntMat$lonAxis)
mecLocs$pos$circleSize <- ceil(10 * mecLocs$pos$coveredAs /
                                      max(mecLocs$pos$coveredAs))
cat(sprintf("uncovered antennas: %d", sum(mecLocs$antennas$MEC == -1)))


# COmpare CDFs
plot(x=cdfDig, pch=19, col="red", xlab=TeX("C_A"), ylab=TeX("$P(X)$"),
     main = "dig-out vs antenna-red CDFs")
lines(x=cdfAntennaRed, pch=18, col="blue", lty=2)
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
plot(x=cdfDig,col="red", legend="dig-out")
lines(x=cdfAntennaRed,col="blue", legend="antenna-reduce")

# Plot last obtained MEC locations
ggmap(map) + 
  ggplot2::stat_contour(data = loggedMecLocsMat,
                        aes(z=intensity, fill=..level..), geom = "polygon",
                        alpha=0.3) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  geom_point(data=mecLocs$pos, aes(x=lon, y=lat, size=circleSize), shape="circle") +
  labs(size = TeX(sprintf("$ceil\\left(\\frac{10·C_A(lon,lat)}{%d} \\right)$",
                              max(mecLocs$pos$coveredAs))))


# Too computational intensive (~ 6h in Madrid-center)
# timestamp()
# mecIntF <- mecGenInt(lons = allAntennas$lon, lats = allAntennas$lat,
#           radios = allAntennas$radio, maxDiss = maxDiss)
# mecIntFMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
#                         longis = c(region$bl$lon, region$tr$lon),
#                         latisLen = LAT_SAMPLES,
#                         longisLen = LON_SAMPLES, intF = mecIntF) 
# timestamp()




# Decide if the logged representation is wanted
loggedMecFrame <- data.frame(lat = mecFrame$lat, lon = mecFrame$lon,
                             intensity = log(1 + mecFrame$intensity))
# OJO: el 3.8 esta puesto como el valor que te quita los bordes
loggedMecFrame$intensity[loggedMecFrame$intensity < 3.8] <- 3.8

ggmap(map) + 
  ggplot2::stat_contour(data=loggedMecFrame, aes(z=intensity, fill=..level..),
               geom = "polygon", alpha=0.3) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  geom_point(data=mecLocs$pos, aes(x=lon, y=lat)) +
  labs(fill = TeX("$log(1 + \\lambda (lon, lat))$")) 
