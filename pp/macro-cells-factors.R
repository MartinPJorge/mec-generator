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
PEOPLE_MAT <- "../data/people/Madrid-centro/people-lambda"
PEOPLE_LONS <- "../data/people/Madrid-centro/people-longitudes"
PEOPLE_LATS <- "../data/people/Madrid-centro/people-latitudes"
OUT_SQUARE_FACTORS <- NULL
CLI <- FALSE # flag to tell if file is executed from CLI

# How big are the cells where we now how many AAUs are
# according to Luca Cominardi, we have 12 AAUs per square kilometer
SQUARE_SIDE <- 1 # expresed in kilometers
SQUARE_AAUs <- 12
REPULSION <- 1 / (2*ceil(sqrt(SQUARE_AAUs))) # expressed in kilometers


# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 5) {
    stop(paste("Arguments to receive are: ",
      "peopleMat longitudeSamp latitudeSamp AAUs outSquareFactors"))
  } else {
    CLI <- TRUE
    PEOPLE_MAT <- args[1]
    PEOPLE_LONS <- args[2]
    PEOPLE_LATS <- args[3]
    SQUARE_AAUs <- as.numeric(args[4])
    OUT_SQUARE_FACTORS <- args[5]
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
peopleIntpF <- genInterpInt(popMat, lonAxis, latAxis)


# Divide regions in areas of ~ 1km x 1km
regionSquares1 <- divideInSquares(latT = latT, latB = latB,
                                  lonL = lonL, lonR = lonR,
                                  area = SQUARE_SIDE * SQUARE_SIDE)

# Get the mapI MatternII average number of points at each square
avgPoints <- c()
for (sq in 1:nrow(regionSquares1)) {
  square <- regionSquares1[sq,]
  avgSq <- matternIImapIintMapprox3(lonL = square$lonL, lonR = square$lonR,
                                    latB = square$latB, latT = square$latT,
                                    r = REPULSION * 1000,
                                    intensityF = peopleIntpF)
  avgPoints <- c(avgPoints, avgSq)
}
regionSquares1 <- data.frame(lonL = regionSquares1$lonL,
                             lonR = regionSquares1$lonR,
                             latB = regionSquares1$latB,
                             latT = regionSquares1$latT,
                             smaller = regionSquares1$smaller,
                             avgAAUs = avgPoints)


# FInd the factor to multiply intensity function at each square
baseSideSquares <- ceil(sqrt(SQUARE_AAUs))
sideSquares <- seq(from = baseSideSquares, to = baseSideSquares + 5, by = 1)
intFactors <- c(SQUARE_AAUs, SQUARE_AAUs * 2, SQUARE_AAUs * 4, SQUARE_AAUs * 8)
squareFactors <- c()
squareSideSqs <- c()
squareAvgs <- c()

for (row in 1:nrow(regionSquares1)) {
  print("Square")
  print(row)
  currAvg <- Inf
  currFactor <- Inf
  currSquares <- Inf
  square <- regionSquares1[row,]
  win <- owin(xrange = c(square$lonL, square$lonR),
              yrange = c(square$latB, square$latT))
  
  
  # Leave the smaller squares without calculation
  if (square$smaller) {
    squareFactors <- c(squareFactors, Inf)
    squareSideSqs <- c(squareSideSqs, Inf)
    squareAvgs <- c(squareAvgs, Inf)
  } else {
    # Find best setup for the square to have avg number of points greater than
    # SQUARE_AAUs
    for (squares in baseSideSquares) {
      repulsion <- 1 / (2*squares) # expressed in kilometers
      print("  #suqres")
      print(squares)
      
      for (factor_ in intFactors) {
        print("  factor")
        print(factor_)
        lambdaSq <- function(lon, lat) {
          return(factor_ / square$avgAAUs * peopleIntpF(lon, lat))
        }
        
        ps <- c()
        for (i in 1:40) {
          points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win,
                                       r = repulsion * 1000)
          ps <- c(ps, points_$n)
        }
        avgPs <- mean(ps)
        
        # Found avg > SQUARE_AAUs
        if(avgPs > SQUARE_AAUs & avgPs < currAvg) {
          currAvg <- avgPs
          currFactor <- factor_
          currSquares <- squares
        }
      }
    }
    
    # Append the best obtained values
    squareFactors <- c(squareFactors, currFactor)
    squareSideSqs <- c(squareSideSqs, currSquares)
    squareAvgs <- c(squareAvgs, currAvg)
    print("found average")
    print(currAvg)
  }
}


# Substitute the smaller squares parameters by the obtained in the previous
# non-smaller
lastCorrectIdx <- 1
for (i in 1:length(squareAvgs)) {
  if (squareAvgs[i] == Inf) {
    squareFactors[i] <- squareFactors[lastCorrectIdx]
    squareSideSqs[i] <- squareSideSqs[lastCorrectIdx]
  } else {
    lastCorrectIdx <- i
  }
}

regionSquares1 <- data.frame(lonL = regionSquares1$lonL,
                             lonR = regionSquares1$lonR,
                             latB = regionSquares1$latB,
                             latT = regionSquares1$latT,
                             smaller = regionSquares1$smaller,
                             avgAAUs = regionSquares1$avgAAUs,
                             factor = squareFactors,
                             sideSquares = squareSideSqs,
                             avgAAUs = squareAvgs)

if (CLI) {
  write.csv(regionSquares1, file = OUT_SQUARE_FACTORS)
}


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