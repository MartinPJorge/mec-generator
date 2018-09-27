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

########### MADRID PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
MEC_GEN_PARAMS <- "../data/mec-gen-params.json"
LAT_SAMPLES <- 100
LON_SAMPLES <- 100
NUM_MECs <- NULL
MACRO_CELLS_CSV <- "../data/antennas/Madrid-center/macro-cells.csv"
MACRO_CELLS_DIR <- "../data/antennas/Madrid-center/antennasDir"
FEMTO_CELLS_CSV <- "../data/antennas/Madrid-center/femto-cells.csv"
MEC_LOCATIONS_CSV <- "../data/antennas/Madrid-center/MEC-locations.csv"
CLI <- FALSE # flag to tell if file is executed from CLI
METHOD <- "basicm1m2"
RADIO_TECH <- "FDD30kHz2sSPS"
INT_MAT_M1 <- NULL
INT_MAT_M2 <- NULL
OPERATOR <- 7
#####################################


########### COBO CALLEJA PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Cobo-Calleja"
MEC_GEN_PARAMS <- "../data/mec-gen-params.json"
LAT_SAMPLES <- 100
LON_SAMPLES <- 100
NUM_MECs <- NULL
MACRO_CELLS_CSV <- paste("../data/antennas/Cobo-Calleja/",
  "antennas-gen-12AAUs-factor13/Cobo-Calleja-12AAUs-factor13-1.csv", sep = "")
FEMTO_CELLS_CSV <- NULL
MEC_LOCATIONS_CSV <- "../data/antennas/Cobo-Calleja/MEC-locations.csv"
CLI <- FALSE # flag to tell if file is executed from CLI
METHOD <- "basicm1m2"
RADIO_TECH <- "FDD30kHz2sSPS"
INT_MAT_M1 <- NULL
INT_MAT_M2 <- NULL
OPERATOR <- 7
###########################################

########### HOCES DEL CABRIEL PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Hoces-del-Cabriel"
MEC_GEN_PARAMS <- "../data/mec-gen-params.json"
LAT_SAMPLES <- 100
LON_SAMPLES <- 300
NUM_MECs <- NULL
MACRO_CELLS_CSV <- NULL
FEMTO_CELLS_CSV <- NULL
ROAD_CELLS_CSV <- "../data/antennas/Hoces-del-Cabriel/road-antennas.csv"
MEC_LOCATIONS_CSV <- paste("../data/mec-pops/Hoces-del-Cabriel/",
                           "mec-locations-i-iii-road/MEC-locations.csv", sep = "")
CLI <- FALSE # flag to tell if file is executed from CLI
METHOD <- "basicm1m2"
RADIO_TECH <- "FDD30kHz2sSPS"
INT_MAT_M1 <- NULL
INT_MAT_M2 <- NULL
OPERATOR <- 7
###############################################



# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 13) {
    stop(paste("Arguments to receive are: ",
      "regionID #latSamples #lonSamples mecParamsJSON numMECs|<0",
      "radioTech[FDD30kHz2sSPS, TDD120kHz7sSPS, FDD120kHz7sSPS]",
      "method[basicm1, basicm2, dig-outm1, dig-outm2, basicm1m2, dig-outm1m2]",
      "macro-cells.csv|NULL femto-cells.csv|NULL road-cells.csv|NULL",
      "MEClocations.csv intMatM1.csv|NULL intMatM2.csv|NULL"))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    LAT_SAMPLES <- as.numeric(args[2])
    LON_SAMPLES <- as.numeric(args[3])
    MEC_GEN_PARAMS <- args[4]
    NUM_MECs <- as.numeric(args[5])
    if (NUM_MECs < 0) {
      NUM_MECs = NULL
    }
    RADIO_TECH <- args[6]
    METHOD <- args[7]
    MACRO_CELLS_CSV <- args[8]
    if (MACRO_CELLS_CSV == "NULL") {
      MACRO_CELLS_CSV = NULL
    }
    FEMTO_CELLS_CSV <- args[9]
    if (FEMTO_CELLS_CSV == "NULL") {
      FEMTO_CELLS_CSV = NULL
    }
    ROAD_CELLS_CSV <- args[10]
    if (ROAD_CELLS_CSV == "NULL") {
      ROAD_CELLS_CSV = NULL
    }
    MEC_LOCATIONS_CSV <- args[11]
    INT_MAT_M1 <- args[12]
    INT_MAT_M2 <- args[13]
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
femtoCells <- NULL
if (!is.null(FEMTO_CELLS_CSV)) {
  femtoCells <- read.csv(file = FEMTO_CELLS_CSV)
}
roadCells <- NULL
if (!is.null(ROAD_CELLS_CSV)) {
  roadCells <- read.csv(file = ROAD_CELLS_CSV)
}
macroCells <- NULL
if (!is.null(MACRO_CELLS_CSV)) {
  macroCells <- read.csv(file = MACRO_CELLS_CSV)
}
lons <- c()
lats <- c()
radio <- c()
net <- c()
# Join femto-cells and macro-cells
if (!is.null(FEMTO_CELLS_CSV)) {
  for (row in 1:nrow(femtoCells)) {
    antenna <- femtoCells[row,]
    lons <- c(lons, antenna$lon)
    lats <- c(lats, antenna$lat)
    radio <- c(radio, as.character(antenna$radio))
    net <- c(net, antenna$net)
  }
}
# Join road cells and macro-cells
if (!is.null(ROAD_CELLS_CSV)) {
  for (row in 1:nrow(roadCells)) {
    antenna <- roadCells[row,]
    lons <- c(lons, antenna$lon)
    lats <- c(lats, antenna$lat)
    radio <- c(radio, as.character(antenna$radio))
    net <- c(net, 7)
  }
}
if (!is.null(MACRO_CELLS_CSV)) {
  for (row in 1:nrow(macroCells)) {
    macroCell <- macroCells[row,]
    lons <- c(lons, macroCell$lon)
    lats <- c(lats, macroCell$lat)
    radio <- c(radio, as.character(macroCell$radio))
    net <- c(net, OPERATOR)
  }
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

# Write the M1 and M2 matrixes
if (CLI) {
  cat("Storing the M1/M2 matrix\n")
  if (METHOD == "basicm1m2" | METHOD == "dig-outm1m2" ) {
    write.csv(mecIntMats$m1$matrix, file = INT_MAT_M1)
    write.csv(mecIntMats$m2$matrix, file = INT_MAT_M2)
  } else if (METHOD == "basicm1" | METHOD == "dig-outm1") {
    write.csv(mecIntMats$m1$matrix, file = INT_MAT_M1)
  } else if (METHOD == "basicm2" | METHOD == "dig-outm2") {
    write.csv(mecIntMats$m2$matrix, file = INT_MAT_M2)
  }
}




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
cat("Storing the MEC locations\n")
write.csv(mecLocs$pos, file = MEC_LOCATIONS_CSV)



#############################################################
################ CODE BEFORE SCRIPT PURPOSES ################
#############################################################
### # Get the CDF of covered antennas
### cdfAntennaRed <- ecdf(mecLocs$pos$coveredAs)
### medAntRed <- median(mecLocs$pos$coveredAs)
### 
### 
### loggedMecLocsMat <- intMat2Frame(intMat = log(1 + mecLocs$modMat$m1),
###                          latAxis = mecIntMats$m1$latAxis,
###                          lonAxis = mecIntMats$m1$lonAxis)
### mecLocs$pos$circleSize <- ceil(10 * mecLocs$pos$coveredAs /
###                                       max(mecLocs$pos$coveredAs))
### cat(sprintf("uncovered antennas: %d", sum(mecLocs$antennas$MEC == -1)))
### 
### # Get the map
### mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
###                right = region$tr$lon, top = region$tr$lat)
### map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
###               maptype = region$plotDetails$mapType)
### 
### # COmpare CDFs
### plot(x=cdfDig, pch=19, col="red", xlab=TeX("C_A"), ylab=TeX("$P(X)$"),
###      main = "dig-out vs antenna-red CDFs")
### lines(x=cdfAntennaRed, pch=18, col="blue", lty=2)
### legend("bottomright", legend=c("dig-out", "red-antenna"),
###        col=c("red", "blue"), lty=1:2, cex=0.8)
### plot(x=cdfDig,col="red", legend="dig-out")
### lines(x=cdfAntennaRed,col="blue", legend="antenna-reduce")
### 
### # Plot last obtained MEC locations
### ggmap(map) + 
###   ggplot2::stat_contour(data = loggedMecLocsMat,
###                         aes(z=intensity, fill=..level..), geom = "polygon",
###                         alpha=0.3) +
###   scale_fill_gradient(low = "gray", high = "black") +
###   scale_colour_gradient(low = "gray", high = "black") +
###   geom_point(data=mecLocs$pos, aes(x=lon, y=lat, size=circleSize), shape="circle") +
###   labs(size = TeX(sprintf("$ceil\\left(\\frac{10·C_A(lon,lat)}{%d} \\right)$",
###                               max(mecLocs$pos$coveredAs))))
### 
### 
### # Too computational intensive (~ 6h in Madrid-center)
### # timestamp()
### # mecIntF <- mecGenInt(lons = all5gAntennas$lon, lats = all5gAntennas$lat,
### #           radios = all5gAntennas$radio, maxDiss = maxDiss)
### # mecIntFMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
### #                         longis = c(region$bl$lon, region$tr$lon),
### #                         latisLen = LAT_SAMPLES,
### #                         longisLen = LON_SAMPLES, intF = mecIntF) 
### # timestamp()
### 
### 
### 
### 
### # Decide if the logged representation is wanted
### loggedMecFrame <- data.frame(lat = mecFrame$lat, lon = mecFrame$lon,
###                              intensity = log(1 + mecFrame$intensity))
### # OJO: el 3.8 esta puesto como el valor que te quita los bordes
### loggedMecFrame$intensity[loggedMecFrame$intensity < 3.8] <- 3.8
### 
### ggmap(map) + 
###   ggplot2::stat_contour(data=loggedMecFrame, aes(z=intensity, fill=..level..),
###                geom = "polygon", alpha=0.3) +
###   scale_fill_gradient(low = "gray", high = "black") +
###   scale_colour_gradient(low = "gray", high = "black") +
###   geom_point(data=mecLocs$pos, aes(x=lon, y=lat)) +
###   labs(fill = TeX("$log(1 + \\lambda (lon, lat))$")) 
### 
### 
### # Plot the MEC PoPs locations
### shapes <- c()
### for (row in 1:nrow(mecLocs$pos)) {
###   shapes <- c(shapes, ifelse(mecLocs$pos$ring[row] == "M1", yes = 18, no = 15))
### }
### mecLocs$pos$shapes <- as.factor(shapes)
### ggmap(map) + 
###   geom_point(data=mecLocs$pos, aes(x=lon, y=lat, shape=shapes,
###                                    size=circleSize)) +
###   labs(size = TeX("$\\frac{10·C_A(mp)}{max_{mp \\in PoPs} C_A(mp)}$")) +
###   scale_shape_manual(name = "Network ring",
###                      labels = c("M2", "M1"),
###                      values = c(15, 18))
### 
### mecLocs$pos$circleSize <- ceil(10 * mecLocs$pos$coveredAs /
###                                       max(mecLocs$pos$coveredAs))
### 
