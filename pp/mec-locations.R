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
NUM_MECs <- NULL
MACRO_CELLS_CSV <- "../data/antennas/Madrid-center/macro-cells.csv"
ANTENNAS_CSV <- "../data/antennas/Madrid-center/TEL-and-femtos.csv" # femto cells and LTE
OPERATOR <- 7 # net id of Telefonica
MEC_LOCATIONS_CSV <- "../data/antennas/Madrid-center/MEC-locations.csv"
CLI <- FALSE # flag to tell if file is executed from CLI
METHOD <- "basicm1m2"
RADIO_TECH <- "FDD30kHz2sSPS"



# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 11) {
    stop(paste("Arguments to receive are: ",
      "regionID latSamples lonSamples mecParamsJSON numMECs",
      "radioTech[FDD30kHz2sSPS, TDD120kHz7sSPS, FDD120kHz7sSPS]",
      "method[basicm1, basicm2, dig-outm1, dig-outm2, basicm1m2, dig-outm1m2]",
      "macro-cells.csv LTEandFEMTO.csv operatorNET MEClocations.csv"))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    MEC_GEN_PARAMS <- args[2]
    LAT_SAMPLES <- args[3]
    LON_SAMPLES <- arg[4]
    NUM_MECs <- as.numeric(args[5])
    NUM_MECs <- ifelse(NUM_MECs < 0, yes = NULL, no = NUM_MECs)
    RADIO_TECH <- args[6]
    METHOD <- args[7]
    MACRO_CELLS_CSV <- args[8]
    ANTENNAS_CSV <- args[9]
    OPERATOR <- as.numeric(args[10])
    MEC_LOCATIONS_CSV <- args[11]
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
  if (antenna$radio != "LTE") {
    lons <- c(lons, antenna$lon)
    lats <- c(lats, antenna$lat)
    radio <- c(radio, as.character(antenna$radio))
    net <- c(net, antenna$net)
  }
}
for (row in 1:nrow(macroCells)) {
  macroCell <- macroCells[row,]
  lons <- c(lons, macroCell$lon)
  lats <- c(lats, macroCell$lat)
  radio <- c(radio, as.character(macroCell$radio))
  net <- c(net, OPERATOR)
}
all5gAntennas <- data.frame(radio = radio, net = net, lon = lons, lat = lats)


##################################
# Obtain the intensity functions #
##################################
mecIntMat <- NULL
mecIntMats <- NULL
maxDiss <- NULL
maxDisss <- NULL
loggedMecInt <- NULL
loggedMecIntM1 <- NULL
loggedMecIntM2 <- NULL

# Get generation parameters depending on the antennas radio technologies
if (RADIO_TECH == "FDD30kHz2sSPS") {
  mecParams <- list(LTE = mecParams$FDD30kHz2sSPS$LTE,
                    macroCell = mecParams$FDD30kHz2sSPS$macroCell,
                    femtoCell = mecParams$FDD30kHz2sSPS$femtoCell,
                    maxAntennas = mecParams$maxAntennas)
} else if (RADIO_TECH == "TDD120kHz7sSPS") {
  mecParams <- list(LTE = mecParams$TDD120kHz7sSPS$LTE,
                    macroCell = mecParams$TDD120kHz7sSPS$macroCell,
                    femtoCell = mecParams$TDD120kHz7sSPS$femtoCell,
                    maxAntennas = mecParams$maxAntennas)
} else if (RADIO_TECH == "FDD120kHz7sSPS") {
  mecParams <- list(LTE = mecParams$FDD120kHz7sSPS$LTE,
                    macroCell = mecParams$FDD120kHz7sSPS$macroCell,
                    femtoCell = mecParams$FDD120kHz7sSPS$femtoCell,
                    maxAntennas = mecParams$maxAntennas)
}

if (METHOD == "basicm1m2" | METHOD == "dig-outm1m2") {
  # Create the maximum distances for both M1 and M2 network rings
  maxDissM1 <- list(LTE = mecParams$LTE$maxDistanceM1,
                  macroCell = mecParams$macroCell$maxDistanceM1,
                  femtoCell = mecParams$femtoCell$maxDistanceM1)
  maxDissM2 <- list(LTE = mecParams$LTE$maxDistanceM2,
                  macroCell = mecParams$macroCell$maxDistanceM2,
                  femtoCell = mecParams$femtoCell$maxDistanceM2)
  maxDisss <- list(m1 = maxDissM1, m2 = maxDissM2)
  
  
  # Generate the intensity functions corresponding to M1 and M2 network rings
  mecIntMatM1 <- mecIntManhattan(lonL = region$bl$lon, lonR = region$tr$lon,
                  latB = region$bl$lat, latT = region$tr$lat,
                  lonSamples = LON_SAMPLES, latSamples = LAT_SAMPLES,
                  antLons = all5gAntennas$lon, antLats = all5gAntennas$lat,
                  antRadios = as.character(all5gAntennas$radio),
                  maxDiss = maxDissM1)
  mecIntMatM2 <- mecIntManhattan(lonL = region$bl$lon, lonR = region$tr$lon,
                  latB = region$bl$lat, latT = region$tr$lat,
                  lonSamples = LON_SAMPLES, latSamples = LAT_SAMPLES,
                  antLons = all5gAntennas$lon, antLats = all5gAntennas$lat,
                  antRadios = as.character(all5gAntennas$radio),
                  maxDiss = maxDissM2)
  mecIntMats <- list(m1 = mecIntMatM1, m2 = mecIntMatM2)
  
  # Log the intensity matrixes
  loggedMecIntM1 <- intMat2Frame(intMat = log(1 + mecIntMatM1$matrix),
                           latAxis = mecIntMatM1$latAxis,
                           lonAxis = mecIntMatM1$lonAxis)
  loggedMecIntM2 <- intMat2Frame(intMat = log(1 + mecIntMatM2$matrix),
                           latAxis = mecIntMatM2$latAxis,
                           lonAxis = mecIntMatM2$lonAxis)
} else {
  # Obtain the M1 or M2 distances antenna-MEC
  if (METHOD == "basicm1" | METHOD == "digoutm1"){
    maxDiss <- list(LTE = mecParams$LTE$maxDistanceM1,
                    macroCell = mecParams$macroCell$maxDistanceM1,
                    femtoCell = mecParams$femtoCell$maxDistanceM1)
  } else {
    maxDiss <- list(LTE = mecParams$LTE$maxDistanceM2,
                    macroCell = mecParams$macroCell$maxDistanceM2,
                    femtoCell = mecParams$femtoCell$maxDistanceM2)
  }
  
  # Generate the intensity functions corresponding to M1 and M2 network rings
  mecIntMat <- mecIntManhattan(lonL = region$bl$lon, lonR = region$tr$lon,
                  latB = region$bl$lat, latT = region$tr$lat,
                  lonSamples = LON_SAMPLES, latSamples = LAT_SAMPLES,
                  antLons = all5gAntennas$lon, antLats = all5gAntennas$lat,
                  antRadios = as.character(all5gAntennas$radio),
                  maxDiss = maxDiss)
  
  # Log the intensity matrixes
  loggedMecInt <- intMat2Frame(intMat = log(1 + mecIntMat$matrix),
                           latAxis = mecIntMat$latAxis,
                           lonAxis = mecIntMat$lonAxis)
}

# Get the map
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
              maptype = region$plotDetails$mapType)

# Plot the intensity function
# TODO - when m1m2 METHOD is selected, only one can be drawn
ggmap(map) + 
  ggplot2::stat_contour(data = loggedMecInt, aes(z=intensity, fill=..level..),
               geom = "polygon", alpha=0.3) +
  scale_fill_gradient(low = "gray", high = "black") +
  scale_colour_gradient(low = "gray", high = "black") +
  labs(fill = TeX("$log(1 + \\lambda (lon, lat))$")) 



###################################
# Obtain the MEC server locations #
###################################
if (METHOD == "basicm1" | METHOD == "basicm2") {
  mecLocs <- mecLocationAntenna(mecIntMatrix = mecIntMat$matrix,
                                lonAxis = mecIntMat$lonAxis,
                                latAxis = mecIntMat$latAxis, maxDiss = maxDiss,
                                numMECs = NUM_MECs, antLons = all5gAntennas$lon,
                                antLats = all5gAntennas$lat,
                                antRadios = as.character(all5gAntennas$radio),
                                lonL = region$bl$lon, lonR = region$tr$lon,
                                latB = region$bl$lat, latT = region$tr$lat) 
} else if (METHOD == "dig-outm1" | METHOD == "dig-outm2") {
  mecLocs <- mecLocationDig(mecIntMatrix = mecIntMat$matrix,
                                lonAxis = mecIntMat$lonAxis,
                                latAxis = mecIntMat$latAxis, maxDiss = maxDiss,
                                numMECs = NUM_MECs, antLons = all5gAntennas$lon,
                                antLats = all5gAntennas$lat,
                                antRadios = as.character(all5gAntennas$radio),
                                lonL = region$bl$lon, lonR = region$tr$lon,
                                latB = region$bl$lat, latT = region$tr$lat) 
} else if (METHOD == "basicm1m2") {
  mats <- list(m1 = mecIntMats$m1$matrix, m2 = mecIntMats$m2$matrix)
  mecLocs <- mecLocationAntennaM1M2(mecIntMatrixs = mats,
                                lonAxis = mecIntMats$m1$lonAxis,
                                latAxis = mecIntMats$m1$latAxis,
                                maxDisss = maxDisss, numMECs = NUM_MECs,
                                antLons = all5gAntennas$lon,
                                antLats = all5gAntennas$lat,
                                antRadios = as.character(all5gAntennas$radio),
                                maxAnts = list(m1 = mecParams$maxAntennas$M1,
                                               m2 = mecParams$maxAntennas$M2),
                                lonL = region$bl$lon, lonR = region$tr$lon,
                                latB = region$bl$lat, latT = region$tr$lat) 
} else if (METHOD == "dig-outm1m2") {
  mecIntMats <- list(m1 = mecIntMatM1, m2 = mecIntMatM2)
  mats <- list(m1 = mecIntMats$m1$matrix, m2 = mecIntMats$m2$matrix)
  mecLocs <- mecLocationDigM1M2(mecIntMatrix = mats,
                                lonAxis = mecIntMats$m1$lonAxis,
                                latAxis = mecIntMats$m1$latAxis,
                                maxDisss = maxDisss, numMECs = NUM_MECs,
                                antLons = all5gAntennas$lon,
                                antLats = all5gAntennas$lat,
                                antRadios = as.character(all5gAntennas$radio),
                                maxAnts = list(m1 = mecParams$maxAntennas$M1,
                                               m2 = mecParams$maxAntennas$M2),
                                lonL = region$bl$lon, lonR = region$tr$lon,
                                latB = region$bl$lat, latT = region$tr$lat) 
}

# Get the CDF of covered antennas
cdfAntennaRed <- ecdf(mecLocs$pos$coveredAs)
medAntRed <- median(mecLocs$pos$coveredAs)


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
legend("bottomright", legend=c("dig-out", "red-antenna"),
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
# mecIntF <- mecGenInt(lons = all5gAntennas$lon, lats = all5gAntennas$lat,
#           radios = all5gAntennas$radio, maxDiss = maxDiss)
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
