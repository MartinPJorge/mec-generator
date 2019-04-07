##### library(pracma)
##### library(SDMTools)
##### library(ggmap)
##### library(spatstat)
##### library(graphics) 
##### library(rjson)
##### library(cluster)
##### library(latex2exp)
##### library(metR)
##### 
##### REGIONS <- "../data/regions.json"
##### REGION_NAME <- "Madrid-centro"
##### PEOPLE_MAT <- "../data/people/Madrid-centro/people-intensity-matrix.csv"
##### PEOPLE_LONS <- "../data/people/Madrid-centro/people-intensity-longitudes"
##### PEOPLE_LATS <- "../data/people/Madrid-centro/people-intensity-latitudes"
##### OUT_SQUARE_FACTORS <- NULL
##### CLI <- FALSE # flag to tell if file is executed from CLI
##### 
##### # How big are the cells where we now how many AAUs are
##### # according to Luca Cominardi, we have 12 AAUs per square kilometer
##### SQUARE_SIDE <- 1 # expresed in kilometers
##### SQUARE_AAUs <- 12
##### REPULSION <- 1 / (2*ceil(sqrt(SQUARE_AAUs))) # expressed in kilometers
##### 
##### 
##### # Parse arguments if existing to change default global variables
##### args <- commandArgs(trailingOnly=TRUE)
##### if(length(args) > 0)
#####   if(length(args) != 5) {
#####     stop(paste("Arguments to receive are: ",
#####       "peopleMat longitudeSamp latitudeSamp AAUs outSquareFactors"))
#####   } else {
#####     CLI <- TRUE
#####     PEOPLE_MAT <- args[1]
#####     PEOPLE_LONS <- args[2]
#####     PEOPLE_LATS <- args[3]
#####     SQUARE_AAUs <- as.numeric(args[4])
#####     OUT_SQUARE_FACTORS <- args[5]
#####   }
##### 
##### 
##### # Read the people lambda function and the longitude and latitude samples
##### popMat <- as.matrix(read.csv(file = PEOPLE_MAT)[,-1])
##### lonAxis <- scan(file = PEOPLE_LONS)
##### latAxis <- scan(file = PEOPLE_LATS)
##### latT <- max(latAxis)
##### latB <- min(latAxis)
##### lonL <- min(lonAxis)
##### lonR <- max(lonAxis)
##### 
##### # Generate the interpolation function
##### peopleIntpF <- genInterpInt(popMat, lonAxis, latAxis)
##### 
##### 
##### # Divide regions in areas of ~ 1km x 1km
##### regionSquares1 <- divideInSquares(latT = latT, latB = latB,
#####                                   lonL = lonL, lonR = lonR,
#####                                   area = SQUARE_SIDE * SQUARE_SIDE)
##### 
##### # Get the mapI MatternII average number of points at each square
##### avgPoints <- c()
##### for (sq in 1:nrow(regionSquares1)) {
#####   square <- regionSquares1[sq,]
#####   avgSq <- matternIImapIintMapprox3(lonL = square$lonL, lonR = square$lonR,
#####                                     latB = square$latB, latT = square$latT,
#####                                     r = REPULSION * 1000,
#####                                     intensityF = peopleIntpF)
#####   avgPoints <- c(avgPoints, avgSq)
##### }
##### regionSquares1 <- data.frame(lonL = regionSquares1$lonL,
#####                              lonR = regionSquares1$lonR,
#####                              latB = regionSquares1$latB,
#####                              latT = regionSquares1$latT,
#####                              smaller = regionSquares1$smaller,
#####                              avgAAUs = avgPoints)
##### 
##### 
##### # FInd the factor to multiply intensity function at each square
##### baseSideSquares <- ceil(sqrt(SQUARE_AAUs))
##### sideSquares <- seq(from = baseSideSquares, to = baseSideSquares + 5, by = 1)
##### intFactors <- c(SQUARE_AAUs, SQUARE_AAUs * 2, SQUARE_AAUs * 4, SQUARE_AAUs * 8)
##### squareFactors <- c()
##### squareSideSqs <- c()
##### squareAvgs <- c()
##### 
##### for (row in 1:nrow(regionSquares1)) {
#####   cat(sprintf("Square: %d\n", row))
#####   currAvg <- Inf
#####   currFactor <- Inf
#####   currSquares <- Inf
#####   square <- regionSquares1[row,]
#####   win <- owin(xrange = c(square$lonL, square$lonR),
#####               yrange = c(square$latB, square$latT))
#####   
#####   
#####   # Leave the smaller squares without calculation
#####   if (square$smaller) {
#####     squareFactors <- c(squareFactors, Inf)
#####     squareSideSqs <- c(squareSideSqs, Inf)
#####     squareAvgs <- c(squareAvgs, Inf)
#####   } else {
#####     # Find best setup for the square to have avg number of points greater than
#####     # SQUARE_AAUs
#####     for (squares in baseSideSquares) {
#####       repulsion <- 1 / (2*squares) # expressed in kilometers
#####       cat(sprintf("  #squares: %d\n", squares))
#####       
#####       for (factor_ in intFactors) {
#####         cat(sprintf("  factor: %d\n", factor_))
#####         lambdaSq <- function(lon, lat) {
#####           return(factor_ / square$avgAAUs * peopleIntpF(lon, lat))
#####         }
#####         
#####         ps <- c()
#####         for (i in 1:40) {
#####           points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win,
#####                                        r = repulsion * 1000)
#####           ps <- c(ps, points_$n)
#####         }
#####         avgPs <- mean(ps)
#####         
#####         # Found avg > SQUARE_AAUs
#####         if(avgPs > SQUARE_AAUs & avgPs < currAvg) {
#####           currAvg <- avgPs
#####           currFactor <- factor_
#####           currSquares <- squares
#####         }
#####       }
#####     }
#####     
#####     # Append the best obtained values
#####     squareFactors <- c(squareFactors, currFactor)
#####     squareSideSqs <- c(squareSideSqs, currSquares)
#####     squareAvgs <- c(squareAvgs, currAvg)
#####     cat(sprintf("found average: %f\n", currAvg))
#####   }
##### }
##### 
##### 
##### # Substitute the smaller squares parameters by the obtained in the previous
##### # non-smaller
##### lastCorrectIdx <- 1
##### for (i in 1:length(squareAvgs)) {
#####   if (squareAvgs[i] == Inf) {
#####     squareFactors[i] <- squareFactors[lastCorrectIdx]
#####     squareSideSqs[i] <- squareSideSqs[lastCorrectIdx]
#####   } else {
#####     lastCorrectIdx <- i
#####   }
##### }
##### 
##### regionSquares1 <- data.frame(lonL = regionSquares1$lonL,
#####                              lonR = regionSquares1$lonR,
#####                              latB = regionSquares1$latB,
#####                              latT = regionSquares1$latT,
#####                              smaller = regionSquares1$smaller,
#####                              avgAAUs = regionSquares1$avgAAUs,
#####                              factor = squareFactors,
#####                              sideSquares = squareSideSqs,
#####                              avgAAUs = squareAvgs)
##### 
##### if (CLI) {
#####   write.csv(regionSquares1, file = OUT_SQUARE_FACTORS)
##### }


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





############### FROM HER ON IS PACKAGED CODE #########################
library(pracma)
library(SDMTools)
library(spatstat)



#' @description MaternII thinning based on Vicenty's distance
#' @param X labeled point pattern object of spatstat (x=lon, y=lat)
#' @return list(lat, lon, marks, n)
vicentyMatIIthin <- function(X, r) {
  survivorsLat <- c()
  survivorsLon <- c()
  survivorsMarks <- c()
  
  for (i in 1:X$n) {
    look <- X[i]
    
    j <- 1
    foundMinor <- FALSE
    while (!foundMinor & j < X$n) {
      if(i != j) {
        dist <- distanceMulti(lat1 = look$y, lon1 = look$x,
                              lat2 = X[j]$y, lon2 = X[j]$x)
        if (dist < r) {
          foundMinor <- X[j]$marks < look$marks
        }
      }
      j <- j + 1
    }
    
    if (!foundMinor) {
      survivorsLat <- c(survivorsLat, look$y)
      survivorsLon <- c(survivorsLon, look$x)
      survivorsMarks <- c(survivorsMarks, look$marks)
    }
  }
  
  return(list(x = survivorsLon, y = survivorsLat, marks = survivorsMarks,
                    n = length(survivorsLat)))
}


#' @description generates a MaternII process using the markI labeling.
#' Marks of points are higher for lower intensity. The process is for points in
#' maps, and it uses Vicenty's distance for the inhibition.
#' @param lambda intensity function receiving (lon, lat)
#' @param win owin rectangle with the limits
#' @param r MaternII inhibition radius
#' @return data.frame(x, y, marks, n)
jorgeMaternIImapI <- function(lambda, win, r) {
  P <- rpoispp(lambda, win = win, nsim = 1)
  if(P$n == 0)
    return(P)
  
  # Generate the marks for the points
  marks <- c()
  for (i in 1:P$n) {
    marks <- c(marks, markI(lon = P[i]$x, lat = P[i]$y, lambda))
  }
  Plab <- rlabel(P, labels = marks)
  
  return(vicentyMatIIthin(Plab, r = r))
}


#' @description Obtains the sides of a square in a map to get the asked area.
#' @param topLat top latitude for the square
#' @param area area of the square given in square kilometers unit
#' @return list(lonSize, latSize)
findSquare <- function(topLat, area) {
  side <- sqrt(area * 1000^2) # from square kilometers to square meters
  bottomLat <- destination(lat = topLat, lon = 0, bearing = 180,
              distance = side)$lat2
  
  lonSizeTop <- destination(lat = topLat, lon = 0, bearing = 90,
                            distance = side)$lon2
  lonSizeBottom <- destination(lat = bottomLat, lon = 0, bearing = 90,
                            distance = side)$lon2
  lonSize <- (lonSizeTop + lonSizeBottom) / 2
  
  return(list(lonSize = lonSize, latSize = topLat - bottomLat))
}


#' @description Creates an intensity measure function with a matrix of the
#' intensity function values along some longitude and latitude axis.
#' @param intMat matrix with (longitude, latitude) values
#' @param lonAxis vector of longitudes in the matrix
#' @param latAxis vector of latitudes in the matrix
#' @return function that receives (lon, lat) values
#' @note if the interpolator function receives coordinates outside the axis
#' limits, it returns 0
genInterpInt <- function(intMat, lonAxis, latAxis) {
  args <- local({
    list(intMat = intMat, lonAxis = lonAxis, latAxis = latAxis)
  })
  
  interpolator <- function(lon, lat) {
    # interp2 function states y=rows and x=cols, hence y=longitude x=latitude
    interpVal <- interp2(x = args$latAxis, y = args$lonAxis, Z = args$intMat,
            xp = lat, yp = lon, method = "linear")
    # replace NA values for coords outside grid by zero
    interpVal <- replace(interpVal, is.na(interpVal), 0)
    
    return(interpVal)
  }
  
  return(interpolator)
}


#' @description It divides the given rectangle in squares of a certain area.
#' @param latT top latitude of the rectangle
#' @param latB bottom latitude of the rectangle
#' @param lonL left longitude of the rectangle
#' @param lonR right longitude of the rectangle
#' @return data.frame(latT, latB, lonL, lonR, smaller) each row is a square
#' @note it goes top-left to bottom-right
#' @note some squares can be rectangles that fill the region
#' @note it uses findSquare() to obtain the area approximation
#' @note the smaller boolean of the data.frame() specifies if the square is a
#' limiting one, and hence it is smaller than a square kilometer
divideInSquares <- function(latT, latB, lonL, lonR, area) {
  currLat <- latT
  currLon <- lonL
  
  # Grid lines
  lonLines <- c()
  latLines <- c()
  
  areaSquare <- findSquare(topLat = latT, area = area)
  
  # Get the latitude lines of the grid
  while (currLat > latB) {
    latLines <- c(latLines, currLat)
    currLat <- currLat - areaSquare$latSize
  }
  latLines <- c(latLines, latB)
  
  # Get the longitude lines of the grid
  while (currLon < lonR) {
    lonLines <- c(lonLines, currLon)
    currLon <- currLon + areaSquare$lonSize
  }
  lonLines <- c(lonLines, lonR)
  
  lonsL <- c()
  lonsR <- c()
  latsT <- c()
  latsB <- c()
  smaller <- c()
  for (ln in 1:(length(lonLines)-1)) {
    for (lt in 1:(length(latLines) - 1)) {
      lonsL <- c(lonsL, lonLines[ln])
      lonsR <- c(lonsR, lonLines[ln + 1])
      
      latsT <- c(latsT, latLines[lt])
      latsB <- c(latsB, latLines[lt + 1])
      
      if (lt == length(latLines) - 1 | ln == length(lonLines) - 1) {
        smaller <- c(smaller, TRUE)
      } else {
        smaller <- c(smaller, FALSE)
      }
    }
  }
  
  return(data.frame(latT = latsT, latB = latsB,
                    lonL = lonsL, lonR = lonsR, smaller = smaller))
}


#' @description It obtains a mark given an intensity function value
#' @param intFval intensity function value
#' @return 1 / intF(lon, lat)
efMarkI <- function(intFval) {
  mark <- 1 / intFval
  
  # Random number to diferentiate between coordinates with same intensity
  # zeros <- abs(log(mark))
  # rand_ <- runif(1, min = 1e-1*(zeros - 2), max = 1e-1*(zeros - 1))
  
  #return(mark + rand_)
  return(mark)
}

#' @description Marks a point as the inverse of its intensity
#' @param lon longitude coord
#' @param lat latitude coord
#' @param intF intensity function to evaluate (lon, lat) intensity
#' @return 1 / intF(lon, lat)
markI <- function(lon, lat, intF) {
  return(efMarkI(intF(lon, lat)))
}


#' @description Approximates the intensity measure of a mattern II hard core
#' process with markI marks within the rectangle of specified limits.
#' @param lonL left longitude coordinate for the rectangle
#' @param lonR right longitude coordinate for the rectangle
#' @param latB bottom latitude coordinate for the rectangle
#' @param latT top coordinate for the rectangle
#' @param r inhibition radius
#' @param intensityF intensity function associated with the point process
#' @return the Q value of the integral2 function
#' @note the probability of x being the minimum at its ball(r) is considered as
#' Prob(N( A:={x': m(x') < m(x)} ) = 0)
matternIImapIintMapprox3 <- function(lonL, lonR, latB, latT, r, intensityF) {
  integrand4 <- function(lon, lat) {
    centMark <- markI(lon, lat, intensityF)
    
    filterBigger <- function(rho, theta) {
      otherLambdas <- intensityF(lon + rho*cos(theta), lat + rho*sin(theta))
      otherMark <- efMarkI(otherLambdas)
      otherMark <- replace(otherMark, otherMark >= centMark, 0)
      indexFunc <- replace(otherMark, otherMark != 0, 1)
      
      return(otherLambdas * indexFunc * rho)
    }
    
    probNotMin <- integral2(filterBigger, xmin = 0, xmax = r,
                            ymin = 0, ymax = 2*pi)$Q
    
    return(exp(-1 * probNotMin) * intensityF(lon, lat))
  }
  
  return(integral2(integrand4, xmin = lonL, xmax = lonR,
                          ymin = latB, ymax = latT)$Q)
}


#' @description Get the average number of points generated at each square with
#' the used intensity function
#' @param regionSquares data.frame with the square regions information
#' @param repulsion kilometers of distance between cell antennas
#' @param interpIntF interpolated intensity function for the cell antennas
#' generation
#' @return data.frame regionSquares extended with the average number of cells
#' per square
getSquareAvgs <- function(regionSquares, repulsion, interpIntF) {
  # Get the mapI MatternII average number of points at each square
  avgPoints <- c()
  for (sq in 1:nrow(regionSquares)) {
    square <- regionSquares[sq,]
    avgSq <- matternIImapIintMapprox3(lonL = square$lonL, lonR = square$lonR,
                                      latB = square$latB, latT = square$latT,
                                      r = repulsion * 1000,
                                      intensityF = interpIntF)
    avgPoints <- c(avgPoints, avgSq)
  }
  
  return(data.frame(lonL = regionSquares$lonL,
                               lonR = regionSquares$lonR,
                               latB = regionSquares$latB,
                               latT = regionSquares$latT,
                               smaller = regionSquares$smaller,
                               avgCells = avgPoints))
}


#' @description Obtains the factor used to multiply the intensity function at
#' each square within a certain region.
#' @param squareCells number of cell antennas to generate at each square
#' @param regionSquares result of divideInSquares()
#' @param interpIntF interpolation of the intensity function
#' @return named list with the intensity factor per each square, booleans to
#' determin if its a limit square, and the average number of cell antennas that
#' each square generates with the specified factor
#' list(squareFactors, squareSideSqs, squareAvgs)
getSquareFactors <- function(squareCells, regionSquares, interpIntF) {
    
  # FInd the factor to multiply intensity function at each square
  baseSideSquares <- ceil(sqrt(squareCells))
  sideSquares <- seq(from = baseSideSquares, to = baseSideSquares + 5, by = 1)
  intFactors <- c(squareCells, squareCells * 2, squareCells * 4, squareCells * 8)
  squareFactors <- c()
  squareSideSqs <- c()
  squareAvgs <- c()
  
  for (row in 1:nrow(regionSquares)) {
    # catsprintf("Square: %d\n", row))
    currAvg <- Inf
    currFactor <- Inf
    currSquares <- Inf
    square <- regionSquares[row,]
    win <- owin(xrange = c(square$lonL, square$lonR),
                yrange = c(square$latB, square$latT))
    
    
    # Leave the smaller squares without calculation
    if (square$smaller) {
      squareFactors <- c(squareFactors, Inf)
      squareSideSqs <- c(squareSideSqs, Inf)
      squareAvgs <- c(squareAvgs, Inf)
    } else {
      # Find best setup for the square to have avg number of points greater than
      # squareCells
      for (squares in baseSideSquares) {
        repulsion <- 1 / (2*squares) # expressed in kilometers
        # catsprintf("  #squares: %d\n", squares))
        
        for (factor_ in intFactors) {
          # catsprintf("  factor: %d\n", factor_))
          lambdaSq <- function(lon, lat) {
            return(factor_ / square$avgCells* interpIntF(lon, lat))
          }
          
          ps <- c()
          for (i in 1:40) {
            points_ <- jorgeMaternIImapI(lambda = lambdaSq, win = win,
                                         r = repulsion * 1000)
            ps <- c(ps, points_$n)
          }
          avgPs <- mean(ps)
          
          # Found avg > squareCells
          if(avgPs > squareCells & avgPs < currAvg) {
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
      # catsprintf("found average: %f\n", currAvg))
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
  
    return(list(squareFactors = squareFactors, squareSideSqs = squareSideSqs,
                squareAvgs = squareAvgs))
}


#' @description Based on the calculus of baseIntensity(), it determines which
#' factor must multiply each square within the treated region to generate the
#' specified average number of cell antennas.
#' @param squareSide side of a region square in kilometers
#' @param squareCells number of cell antennas within a region square
#' @param intenMatrix intensity matrix returned by baseIntensity()
#' @param lonAxis longitude axis returned by baseIntensity()  
#' @param latAxis latitude axis returned by baseIntensity()
#' @return data.frame with region square longitude and latitude limits,
#' if they are small squares, their intensity factor, number of side squares,
#' and average number of cells that are generated, i.e.,
#' data.frame(lonL, lonR, latB, latT, smaller, factor, sideSquares, avgCells)
intensityFactors <- function(squareSide = 1, squareCells = 12,
                           intenMatrix = NULL, lonAxis = NULL, latAxis = NULL) {
  repulsion <- 1 / (2*ceil(sqrt(squareCells))) # expressed in kilometers
  
  
  latT <- max(latAxis)
  latB <- min(latAxis)
  lonL <- min(lonAxis)
  lonR <- max(lonAxis)
  
  # Generate the interpolation intensity function
  interpIntF <- genInterpInt(intenMatrix, lonAxis, latAxis)
  # Divide regions in areas of ~ (squareSide km x squareSide km), and avg cells
  regionSquares <- divideInSquares(latT = latT, latB = latB,
                                    lonL = lonL, lonR = lonR,
                                    area = squareSide * squareSide)
  regionSquares <- getSquareAvgs(regionSquares, repulsion, interpIntF)
  
  # Get the intensity factors per each square to get the desired number of cells
  squareFacts <- getSquareFactors(squareCells = squareCells,
                                  regionSquares = regionSquares,
                                  interpIntF = interpIntF)
  
  return(data.frame(lonL = regionSquares$lonL,
                    lonR = regionSquares$lonR,
                    latB = regionSquares$latB,
                    latT = regionSquares$latT,
                    smaller = regionSquares$smaller,
                    avgCells = regionSquares$avgCells,
                    factor = squareFacts$squareFactors,
                    sideSquares = squareFacts$squareSideSqs,
                    avgCells = squareFacts$squareAvgs))
}


