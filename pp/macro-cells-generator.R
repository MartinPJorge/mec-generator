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
SQUARE_FACTORS <- "../data/people/Madrid-center/people-lambda"
MACRO_CELLS_CSV <- "../data/antennas/Madrid-center/macro-cells.csv"
PEOPLE_MAT <- "../data/people/Madrid-center/people-lambda"
PEOPLE_LONS <- "../data/people/Madrid-center/people-longitudes"
PEOPLE_LATS <- "../data/people/Madrid-center/people-latitudes"
CLI <- FALSE # flag to tell if file is executed from CLI

# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 5) {
    stop(paste("Arguments to receive are: ",
      "peopleMat longitudeSamp latitudeSamp squareFactors macroCellsCSV"))
  } else {
    CLI <- TRUE
    PEOPLE_MAT <- args[1]
    PEOPLE_LONS <- args[2]
    PEOPLE_LATS <- args[3]
    SQUARE_FACTORS <- args[4]
    MACRO_CELLS_CSV <- args[5]
  }


# Read the people lambda function and the longitude and latitude samples
popMat <- as.matrix(read.csv(file = PEOPLE_MAT)[,-1])
lonAxis <- scan(file = PEOPLE_LONS)
latAxis <- scan(file = PEOPLE_LATS)
latT <- max(latAxis)
latB <- min(latAxis)
lonL <- min(lonAxis)
lonR <- max(lonAxis)

# Generate the interpolation function
print("Generating the people intensity function as input matrix interpolation")
peopleIntpF <- genInterpInt(popMat, lonAxis, latAxis)


squareFactors <- read.csv(file = SQUARE_FACTORS)


# Iterate through each square and generate the antennas
macroLons <- c()
macroLats <- c()
lonLs <- c()
lonRs <- c()
latBs <- c()
latTs <- c()
squareNos <- c()
for (row in 1:nrow(squareFactors)) {
  square <- squareFactors[row,]
  
  cat(sprintf("Generating macro cells for square %d\n", row))
  cat(sprintf("  lonL=%f lonR=%f latB=%f latT=%f\n", square$lonL, square$lonR,
          square$latB, square$latT))
  
  squareWin <- owin(xrange = c(square$lonL, square$lonR),
                    yrange = c(square$latB, square$latT))
  squareRepulsion <- 1 / (2 * square$sideSquares)
  
  lambdaSq <- function(lon, lat) {
    return(square$factor / square$avgAAUs * peopleIntpF(lon, lat))
  }
  
  # Generate the macro-cell antennas
  squareMacroCells <- jorgeMaternIImapI(lambda = lambdaSq, win = squareWin,
                                    r = squareRepulsion * 1000)
  
  cat(sprintf("  %d macro-cells generated\n", squareMacroCells$n))
  
  # Append the antennas and square information
  macroLons <- c(macroLons, squareMacroCells$x)
  macroLats <- c(macroLats, squareMacroCells$y)
  lonLs <- c(lonLs, rep(square$lonL, squareMacroCells$n))
  lonRs <- c(lonRs, rep(square$lonR, squareMacroCells$n))
  latBs <- c(latBs, rep(square$latB, squareMacroCells$n))
  latTs <- c(latTs, rep(square$latT, squareMacroCells$n))
  squareNos <- c(squareNos, rep(row, squareMacroCells$n))
}


# Store the generated antennas
write.csv(data.frame(lon = macroLons, lat = macroLats, squareLonL = lonLs,
                     squareLonR = lonRs, squareLatB = latBs, squareLatT = latTs,
                     squareNum = squareNos,
                     radio = rep("macro-cell", length(lonRs))),
          file = MACRO_CELLS_CSV)


## # TEST FOR A CERTAIN SQUARE
## REPULSION <- 1 / (2*(1 + ceil(sqrt(SQUARE_AAUs)))) # expressed in kilometers
## sqI <- 211
## square <- regionSquares1[sqI,]
## win <- owin(xrange = c(square$lonL, square$lonR),
##             yrange = c(square$latB, square$latT))
## lambdaSq <- function(lon, lat) {
##   # return(12/square$avgAAUs * peopleIntpF(lon, lat))
##   return(24/square$avgAAUs * peopleIntpF(lon, lat))
## }
## timestamp()
## points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win, r = REPULSION * 1000)
## print(points_$n)
## plot(x = points_$x, y = points_$y)
## timestamp()
## 
## ps <- c()
## for (i in 1:40) {
##   print(i)
##   points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win, r = REPULSION * 1000)
##   ps <- c(ps, points_$n)
## }
## print("real average")
## print(mean(ps))
## print("estimation")
## avgSq <- matternIImapIintMapprox3_V(lonL = square$lonL, lonR = square$lonR,
##                                   latB = square$latB, latT = square$latT,
##                                   r = REPULSION * 1000,
##                                   intensityF = lambdaSq)
## print(avgSq)
## 
## 
## ############ TEST FOR EVERY SQUARE
## for (row in 1:nrow(regionSquares1)) {
##   square <- regionSquares1[row,]
##   win <- owin(xrange = c(square$lonL, square$lonR),
##               yrange = c(square$latB, square$latT))
##   lambdaSq <- function(lon, lat) {
##     # return(12/square$avgAAUs * peopleIntpF(lon, lat))
##     return(12/square$avgAAUs * peopleIntpF(lon, lat))
##   }
##   
##   # Get the average
##   
##   ps <- c()
##   for (i in 1:40) {
##     points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win, r = REPULSION * 1000)
##     ps <- c(ps, points_$n)
##   }
##   print("square:")
##   print(row)
##   print("Average AAUs:")
##   print(mean(ps))
##   
##   avgSq <- matternIImapIintMapprox3(lonL = square$lonL, lonR = square$lonR,
##                                     latB = square$latB, latT = square$latT,
##                                     r = REPULSION * 1000,
##                                     intensityF = lambdaSq)
##   print(avgSq)
## }
## ##########################################
## 
## 
## ############ PLOT EVOLUTION AS SCALE INCreASES
## REPULSION <- 1 / (2*(100 + ceil(sqrt(SQUARE_AAUs)))) # expressed in kilometers
## scl <- seq(from=12, by=20, to=300)
## avvgss_ <- c()
## real_vals <- c()
## for (escalado in scl) {
##   print(escalado)
##   lambdaSq <- function(lon, lat) {
##     return(escalado/square$avgAAUs * peopleIntpF(lon, lat))
##   }
##   avvgss_ <- c(avvgss_, matternIImapIintMapprox3(lonL = square$lonL, lonR = square$lonR,
##                                     latB = square$latB, latT = square$latT,
##                                     r = REPULSION * 1000,
##                                     intensityF = lambdaSq))
##   points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win, r = REPULSION * 1000)
##   real_vals <-c(real_vals, points_$n)
## }
## plot(x = scl, y = avvgss_)
## lines(x = scl, y = real_vals)
## ###########################################
## 
## 
## ns <- c()
## for (i in 1:40) {
##   print(i)
##   points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win, r = REPULSION * 1000)
##   ns <- c(ns, points_$n)
## }
## print(mean(ns))
## 
## for (i in 1:points_$n) {
##   for (j in 1:points_$n) {
##     if(i != j) {
##       print(distance(lon1 = points_$x[i], lat1 = points_$y[i],
##                      lon2 = points_$x[j], lat2 = points_$y[j])$distance)
##     }
##   }
## }
## plot(x = points_$x, y = points_$y)
## 
## 
## filled.contour(x = lonAxis, y = latAxis, z = popMat)
## 
## 