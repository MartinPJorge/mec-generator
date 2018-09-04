library(pracma)
library('SDMTools')

#' @description Wrapper for the Vicenty's distance of the SDMTools package. It
#' allows receiving a vector for lat1 and lon1.
#' @param lat1 latitude coordinate/s
#' @param lon1 longitude coordinate/s
#' @param lat2 latitude coordinate
#' @param lon2 longitude coordinate
#' @return distance/s between (lat1,lon1) and (lat2,lon2)
distanceMulti <- function(lat1, lon1, lat2, lon2) {
  distances <- rep(0, length(lat1))
  
  for (i in 1:length(distances)) {
    distances[i] <- SDMTools::distance(lat1 = lat1[i], lon1 = lon1[i],
                                       lat2 = lat2, lon2 = lon2)$distance
  }
  
  return(distances)
}


#' Obtains the coordinates from string values stored in regions for JSONs
#' @param coordsStr "latitude,longitude"
#' @return data.frame(lat, lon)
getCoords <- function(coordsStr) {
  coords <- strsplit(coordsStr, ",")[[1]]
  return(data.frame(lon = as.numeric(coords[2]), lat = as.numeric(coords[1])))
}

#' It generates an intensity function that is the sum of more intensity
#' functions centered at certain locations.
#' Note: the distance used is geodetic WGS84 using Vicenty's formula
#'       hence centers are WGS84 coordinates (latitude, longitude)
#'
#' @param centers dataframe with lon and lat coordinates of centers
#' @param centers_intensity probability density function to place at each center
#' @param factor multiplying factor in the intensity function
#' @param ... arguments for the probability density function
#'
#' @return intensity function that is the sum of the centered ones.
#'         function(latitude, longitude)
#' @note the centers_intensity must be a radial function
gen_manta <- function(centers, centers_intensity, factor, ...) {
  # Defines the intensity functions for each center
  center_ints <- list()
  for (i in 1:dim(centers)[1]) {
    center_ints[[i]] <- local({
      origin <- centers[i,]
      
      function(longi, lati) {
        # Vicenty's distance
        r <- distanceMulti(lat1=lati, lon1=longi,
                 lat2=origin$lat, lon2=origin$lon)
         
        return(centers_intensity(r, ...))
      }
      
    })
  }
  
  # Creates the function that is sum of intensities
  sum_int <- function(longi, lati) {
    sum_prob <- 0
    i <- 0
    for (center_int in center_ints) {
      sum_prob <- sum_prob + center_int(longi, lati)
      i <- i + 1
    }
    return(sum_prob * factor)
  }
  
  return(sum_int)
}


#' Generates a dataframe evaluating an intensity measure function in a grid of
#' coordinates.
#'
#' @param latis vector (latitude_min, latitude_max)
#' @param longis vector (longitude_min, longitude_max)
#' @param latisLen sampling length in latitude axis
#' @param longisLen sampling length in longitude axis
#' @param intF intensity measure function used to generate the values
#'
#' @return a dataframe having latitude, longitude and intensity
genIntFrame <- function(latis, longis, latisLen, longisLen, intF) {
  latiAxes <- seq(from = latis[1], to = latis[2], length = latisLen)
  longiAxes <- seq(from = longis[1], to = longis[2], length = longisLen)
  
  intFrame <- expand.grid(lon = longiAxes, lat = latiAxes)
  intFrame <- data.frame(lon = intFrame$lon, lat = intFrame$lat,
                          intensity = array(dim = nrow(intFrame)))
  for (row in 1:nrow(intFrame)) {
    intFrame$intensity[row] <- intF(intFrame$lon[row], intFrame$lat[row])
  }
  
  return(intFrame)
}


#' Generates a matrix evaluating an intensity measure function in a grid of
#' coordinates. Rows are for latitude and columns for longitudes.
#'
#' @param latis vector (latitude_min, latitude_max)
#' @param longis vector (longitude_min, longitude_max)
#' @param latisLen sampling length in latitude axis
#' @param longisLen sampling length in longitude axis
#' @param intF intensity measure function used to generate the values
#'
#' @return a list with $matrix $latAxis $lonAxis
#' @note the genrated matrix is indexed as (longitude, latitude)
genIntMatrix <- function(latis, longis, latisLen, longisLen, intF) {
  latiAxes <- seq(from = latis[1], to = latis[2], length = latisLen)
  longiAxes <- seq(from = longis[1], to = longis[2], length = longisLen)
  
  intMatrix <- matrix(nrow = longisLen, ncol = latisLen)
  
  for (lat in 1:latisLen) {
    for (lon in 1:longisLen) {
      intMatrix[lon, lat] <- intF(longiAxes[lon], latiAxes[lat])
    }
  }
  
  return(list(
    matrix = intMatrix,
    latAxis = latiAxes,
    lonAxis = longiAxes
  ))
}


#' @description It transforms the output of genIntMatrix() into a data.frame
#' @param intMat intensity matrix obtained with genIntMatrix()
#' @param latAxis latitude axis obtained with genIntMatrix()
#' @param lonAxis longitude axis obtained with genIntMatrix()
#' @return data.frame(lon, lat, intensity)
intMat2Frame <- function(intMat, latAxis, lonAxis) {
  intensities <- c()
  longitudes <- c()
  latitudes <- c()
  
  for (row in 1:nrow(intMat)) {
    intensities <- c(intensities, intMat[row,])
    longitudes <- c(longitudes, rep(lonAxis[row], length(intMat[row,])))
    latitudes <- c(latitudes, latAxis)
  }
  
  return(data.frame(lon = longitudes, lat = latitudes, intensity = intensities))
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


#' @description Matern I thinning. If a point lies inside the ball B((x,y), r),
#' then the probability of removal is 1, otherwise is 0. Hard-core thinning.
#' @param lon longitude coordinate
#' @param lat latitude coordinate
#' @param ...
#'    - r: radius of the ball
#' @return removal probability
rad_thin <- function(lon, lat, ...) {
  arguments <- list(...)
  r <- arguments[1][[1]]
  print("R vvv")
  print(r)
  print(length(lat))
  print(length(lon))
  deletions <- rep(0, length(lon))
  
  print("asa")
  # Loop through every point
  for (i in seq(1, length(lon))) {
    j <- 1
    found_in_rad <- FALSE
    
    # Checks if there is another point closer than 'r' to (lon[i],y[i])
    while(!found_in_rad & j <= length(lon)) {
      dist <- SDMTools::distance(lat1 = lat[i], lon1 = lon[i],
                                 lat2 = lat[j], lon2 = lon[j])$distance
      if (dist != 0 & dist <= r) {
        found_in_rad <- TRUE
      }
      j <- j + 1
    }
    
    deletions[i] <- if(found_in_rad) 0 else 1
  }
  
  return(deletions)
}


#' @description Step function going from zero to 1 between a and b, it is 1
#' between b and c, finally it goes from 1 to zero between c and d. Outside
#' (a,d) its result is 0. Its purpose is to be used as an intensity function.
#' @param a start of the step function
#' @param b start of the 1 value
#' @param c end of the 1 value
#' @param d end of the step function
#' @return the step function evaluated at x
#' @note the function is based on the smoothstep:
#' https://en.wikipedia.org/wiki/Smoothstep
stepIntensity <- function(x, a, b, c, d) {
  smooth3 <- function(x) {
    return(-20*x^7 + 70*x^6 - 84*x^5 + 35*x^4)
  }
   
  stepInt <- function(x, a, b, c, d) {
    if (a <= x & x <= b) {
      return(smooth3((x - a) / (b - a)))
    }
    else if (b <= x & x <= c) {
      return(1)
    }
    else if (c <= x & x <= d) {
      return(smooth3((x - d) / (c - d)))
    }
    else {
      return(0)
    }
  }
  
  return(sapply(X = x, FUN = stepInt, a = a, b = b, c = c, d = d))
}


#' @description Obtains the intensity measure of an inhomogeneous PPP inside a
#' rectangle in longitude and latitude coordinates.
#' @param lon1 left longitude
#' @param lon2 right longitude
#' @param lat1 down latitude
#' @param lat2 up latitudell
#' @param intensityF intensity function receiving (lon, lat)
#' @return the average number of points in the rectangle
inhIntM <- function(lon1, lon2, lat1, lat2, intensityF) {
  return(integral2(intensityF, xmin = lon1, xmax = lon2,
                   ymin = lat1, ymax = lat2)$Q)
}


#' @description Calculates the intensity measure of a function within a ball
#' @param lon ball center longitude coordinate 
#' @param lat ball center latitude coordinate 
#' @param r ball radius
#' @param intensityF intensity function receiving (lon, lat)
#' @return the average number of points in the ball
ballIntM <- function(lon, lat, r, intensityF) {
  shifted <- function(rho, theta) {
    return(rho * intensityF(lon + rho*cos(theta), lat + rho*sin(theta)))
  }
  
  return(integral2(shifted, xmin = 0, xmax = r, ymin = 0, ymax = 2*pi)$Q)
}


#' @description Obtains the distance matrix of a point pattern using Vicenty's
#' distance
#' @param X point pattern of spatstat package
#' @return a matrix of distances
pairVicentyDist <- function(X) {
  distances <- matrix(ncol = X$n, nrow = X$n)
  for (row in 1:nrow(distances)) {
    distances[row, ] <- distanceMulti(X$lon, X$lat, X[row]$lon, X[row]$lat)
  }
  return(distances)
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


#' @description Generates a MaternII process using markI as label method.
#' it is based on the original spatstat rMatternII function
#' @note this version goes slower than jorgeMaternIImapI
rMaternIImapI <- function(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE,
                          ..., nsim=1, drop=TRUE) {
  
  bigbox <- if(stationary) grow.rectangle(win, r) else win
  Y <- rpoispp(kappa, win = bigbox, ..., nsim = nsim, drop=drop)
  age <- markI(Y$x, Y$y, kappa)
  d2 <- pairVicentyDist(Y)
  close <- (d2 <= r)
  ## random order 1:n
  older <- outer(age, age, ">")
  conflict <- close & older
  ## delete <- apply(conflict, 1, any)
  delete <- spatstat.utils::matrowany(conflict)
  Y <- Y[!delete]
  
  if (stationary) Y <- Y[win]
  
  return(Y)
}


#' @description Centered version of the step intensity function stepIntensity()
#' @param x radius variable
#' @param b step width
#' @param c center section width
#' @return number
centStepIntensity <- function(x, b, c) {
  return(stepIntensity(x = x + b + c/2, a = 0, b = b, c = b + c, d = 2*b + c))
}


#' @description Intensity function for the generation of people in a town. It is
#' based in centStepIntensity function.
#' @param lon longitude coordinate of evaluation point
#' @param lat latitude coordinate of evaluation point
#' @param centLon longitude coordinate of the town center
#' @param centLat latitude coordinate of the town center
#' @param b b parameter of the centStepIntensity() function
#' @param c c parameter of the centStepIntensity() function
townIntF <- function(lon, lat, centLon, centLat, b, c) {
  r <- distanceMulti(lat1 = lat, lon1 = lon, lat2 = centLat, lon2 = centLon)
  
  return(centStepIntensity(x = r, b = b, c = c))
}


#' @description Obtains the rectangle having the location circle inscribed
#' @param centLon longitude coordinate of the center circle
#' @param centLat latitude coordinate of the center circle
#' @param radius circle radius in meters
#' @return list(latB, latT, lonL, lonR)
outerRect <- function(centLon, centLat, radius) {
  latT <- destination(lat = centLat, lon = centLon, bearing = 0,
                        distance = radius)$lat2
  latB <- destination(lat = centLat, lon = centLon, bearing = 180,
                        distance = radius)$lat2
  lonR <- destination(lat = centLat, lon = centLon, bearing = 90,
                        distance = radius)$lon2
  lonL <- destination(lat = centLat, lon = centLon, bearing = 270,
                        distance = radius)$lon2
  
  return(list(latB = latB, latT = latT, lonL = lonL, lonR = lonR))
}


#' @description Obtains the intensity measure of a town whose people points is
#' generated using centStepIntensity()
#' @param centLon longitude coordinate for the center of town
#' @param centLat latitude coordinate for the center of town
#' @param b step width
#' @param c center section width
#' @return the intensity measure of the town
townIntM <- function(centLon, centLat, b, c) {
  outRect <- outerRect(centLon, centLat, radius = b + c/2)
  avgPeople <- integral2(townIntF,
                         xmin = outRect$lonL, xmax = outRect$lonR,
                         ymin = outRect$latB, ymax = outRect$latT,
                         centLon = centLon, centLat = centLat,
                         b = b, c = c)$Q
  return(avgPeople)
}


#' @description Generates an intensity function to generate people in towns.
#' It uses centStepIntensity() function to create the intensity.
#' @param townCents list(lat = c(), lon = c()) coordinates with the town center
#' @param townRads radius covering the town
#' @param townPopus vector with the population of each town
#' @param townDisps dispersion factor of people outside towns
#' @note if townDisp[1]=0.5 it means that outside the town, there are
#' 0.5*townRads[1] kilometers more where people can appear.
genPopuManta <- function(townCents, townRads, townPopus, townDisps) {
  townIntensities <- list()
  intMeasures <- c()
  
  # Get the intensity measure for each center
  for (t in 1:length(townRads)) {
    intM <- townIntM(centLon = townCents$lon[t],
                               centLat = townCents$lat[t],
                               b = townRads[t] * townDisps[t], c = townRads[t])
    intMeasures <- c(intMeasures, intM)
  }
  
  # Generate the town intensity functions
  for (t in 1:length(townRads)) {
    townIntensities[[t]] <- local({
      centLat <- townCents$lat[t]
      centLon <- townCents$lon[t]
      rad <- townRads[t]
      disp <- townDisps[t]
      popu <- townPopus[t]
      intMe <- intMeasures[t]
      
      function(lon, lat) {
        r_ <- distanceMulti(lon1 = lon, lat1 = lat,
                            lon2 = centLon, lat2 = centLat)
        b_ <- rad * disp
        c_ <- rad
        return(popu / intMe * centStepIntensity(r_, b_, c_))
      }
    })
  }
  
  # Final function adding all the intensity functions for the towns
  sumTownInt <- function(lon, lat) {
    totalInt <- 0 
    for (townIntensity in townIntensities) {
      totalInt <- totalInt + townIntensity(lon, lat)
    }
    return(totalInt)
  }
  
  return(sumTownInt)
}


#' @description Store the population point coordinates in a CSV
#' @param popuPoints points object from spatstat (x=lon,y=lat)
#' @param outCSV output CSV file where point coordinates are stored
popuToCSV <- function(popuPoints, outCSV) {
  popuCoords <- data.frame(lon = popuPoints$x, lat = popuPoints$y)
  write.csv(popuCoords, file = outCSV)
}


#' @description COST-Hata model for the path loss of a signal
#' @param fc carrier frequency in MHz
#' @param hb BS's height in meters
#' @param hm MS's height in meters
#' @param distance kilometer distance between user and antenna
#' @param metropolitan boolean indicating if it is a dense metropolitan area
#' @return path loss in dBm
COST231Hata <- function(fc, hb, hm, distance, metropolitan = TRUE) {
  a <- function(hm, fc, metropolitan) {
    a_ <- 0
    
    if(!metropolitan)
      a_ <- (1.1 * log10(fc) - 0.7) * hm - 1.56 * log10(fc) - 0.8
    else if(fc < 400)
      a_ <- 8.28 * (log10(1.54 * hm))^2 - 1.1 
    else
      3.2 * (log10(11.75 * hm))^2 - 4.97
    
    return(a_)
  }
  
  A <- 46.3 + 33.9 * log10(fc) - 13.82 * log10(hb) - a(hm)
  B <- 44.9 - 6.55 * log10(hb)
  C <- if(metropolitan) 3 else 0
  
  return(A + B * log10(distance) + C)
}


#' @description COST 231 Walfisch-Ikegami LOS model for path loss of a signal.
#' @param fc carrier frequency in MHz
#' @param distance distance between user and antenna in meters
#' @return path loss in dB
#' @note Appropiate for distances between 0.02-5 km and carrier frequencies
#' between 800 and 2000 MHz. LOS means Line Of Sight, and it assumes no
#' building between receiver and antenna.
COST231WalfischIkegamiLOS <- function(fc, distance) {
  # Distance in formula is in kilometers
  d <- if(distance < 20) 0.02 else distance / 1000
  
  return(42.6 + 26 * log(d) + 20 * log(fc))
}


############## GRID UTILITIES ##############


#' @description Obtains the limits of the i-th division of the (B - A) interval.
#' @param A low interval bound
#' @param B upper interval bound
#' @param divisionL number of divisions of the (B - A) interval
#' @return list(a = low limit, b = upper limit)
getLimits <- function(A, B, divisionL, i) {
  if(B < A | divisionL < 1 | i < 1)
    stop("Inside getLimits() I've received B < A, divisionL < 1 or i < 1")
  return(list(
    a = A + (B - A) / divisionL * (i - 1),
    b = A + (B - A) / divisionL * i
  ))
}


#' @description It obtains the cell inside the grid corresponding to the given
#' longitude and latitude coordinates.
#' @param lon longitude coordinate to locate
#' @param lat latitude coordinate to locate
#' @param lonL longitude left limit
#' @param lonR longitude right limit
#' @param latB latitude botom limit
#' @param latT latitude top limit
#' @param lonLen division length of the longitude grid
#' @param latLen division length of the latitude grid
#' @return list(x, y) numbering of the associated grid
#' @note the grids are numbered bottom up and left to right:
#'  ___ ___ ___
#' | 4 | 5 | 6 |
#' |---|---|---|
#' | 1 | 2 | 3 | < y
#'  --- --- ---
#'   ^
#'   x
getGrid <- function(lon, lat, lonL, lonR, latB, latT, lonLen, latLen) {
  if(lon < lonL | lon > lonR)
    stop("lon coordinate is outside the bounds")
  if(lat < latB | lat > latT)
    stop("lat coordinate is outside the bounds")
  
  lonSepLen <- (lonR - lonL) / lonLen
  lonSquare <- if (lon == lonR) lonLen else floor((lon - lonL) / lonSepLen) + 1
  latSepLen <- (latT - latB) / latLen
  latSquare <- if (lat == latT) latLen else floor((lat - latB) / latSepLen) + 1
  
  return(list(x = lonSquare, y = latSquare))
}

######### TODO - remove, it's substituted by gridData
# #' @description It divides the antennas' data in a grid
# #' @param antennasCSV CSV file where the antennas' data is present
# #' @param region list(bl, br, tl, tr, repulsionRadius, plotDetails, populations)
# #' @param lonN number of divisions in the longitude axis
# #' @param latN number of divisions in the latitude axis
# #' @param griddedJSON output file name where to store the gridded info.
# gridAntennas <- function(antennasCSV, region, lonN, latN, griddedJSON) {
#   antennas <- read.csv(antennasCSV)
#   gridObj <- list(latB = region$bl$lat, latT = region$tr$lat,
#                   lonL = region$bl$lon, lonR = region$tr$lon,
#                   lonN = lonN, latN = latN)
#   
#   # Initialize the squares structure
#   squares <- list()
#   for (loni in 1:lonN) {
#     lonLims <- getLimits(A = gridObj$lonL, B = gridObj$lonR,
#                          divisionL = lonN, i = loni)
#     squares[[loni]] <- list()
#     for (lati in 1:latN) {
#       latLims <- getLimits(A = gridObj$latB, B = gridObj$latT,
#                            divisionL = latN, i = lati)
#       squares[[loni]][[lati]] <- list(latB = latLims$a, latT = latLims$b,
#                                   lonL = lonLims$a, lonR = lonLims$b,
#                                   antennas = list())
#     }
#   }
#   
#   # Put each antenna inside the corresponding square of the grid
#   for (row in 1:nrow(antennas)) {
#     assocGrid <- getGrid(lon = antennas$lon[row], lat = antennas$lat[row],
#                          lonL = gridObj$lonL, lonR = gridObj$lonR,
#                          latB = gridObj$latB, latT = gridObj$latT,
#                          lonLen = gridObj$lonN, latLen = gridObj$latN)
#     squareX <- assocGrid$x
#     squareY <- assocGrid$y
#     numAntennas <- length(squares[[squareX]][[squareY]]$antennas)
#     squares[[squareX]][[squareY]]$antennas[[numAntennas + 1]] <- list(
#       lon = antennas$lon[row], lat = antennas$lat[row],
#       radio = antennas$radio[row]
#     )
#   }
#   
#   # Write info. to a JSON file
#   gridObj$squares <- squares
#   write(toJSON(gridObj, indent = 4, method = "C"), griddedJSON)
# }


#' @description It divides the geographical data in a grid
#' @param inCSV CSV file where the geographical data is present
#' @param region list(bl, br, tl, tr, repulsionRadius, plotDetails, populations)
#' @param lonN number of divisions in the longitude axis
#' @param latN number of divisions in the latitude axis
#' @param griddedJSON output file name where to store the gridded info.
gridData <- function(inCSV, antennas = TRUE, region, lonN, latN, griddedJSON) {
  inData <- read.csv(inCSV)
  gridObj <- list(latB = region$bl$lat, latT = region$tr$lat,
                  lonL = region$bl$lon, lonR = region$tr$lon,
                  lonN = lonN, latN = latN, assignedDis = FALSE)
  
  # Initialize the squares structure
  squares <- list()
  for (loni in 1:lonN) {
    lonLims <- getLimits(A = gridObj$lonL, B = gridObj$lonR,
                         divisionL = lonN, i = loni)
    squares[[loni]] <- list()
    for (lati in 1:latN) {
      latLims <- getLimits(A = gridObj$latB, B = gridObj$latT,
                           divisionL = latN, i = lati)
      if (antennas)
        squares[[loni]][[lati]] <- list(latB = latLims$a, latT = latLims$b,
                                    lonL = lonLims$a, lonR = lonLims$b,
                                    antennas = list())
      else
        squares[[loni]][[lati]] <- list(latB = latLims$a, latT = latLims$b,
                                    lonL = lonLims$a, lonR = lonLims$b,
                                    people = list())
    }
  }
  
  # Put each person inside the corresponding square of the grid
  for (row in 1:nrow(inData)) {
    assocGrid <- getGrid(lon = inData$lon[row], lat = inData$lat[row],
                         lonL = gridObj$lonL, lonR = gridObj$lonR,
                         latB = gridObj$latB, latT = gridObj$latT,
                         lonLen = gridObj$lonN, latLen = gridObj$latN)
    squareX <- assocGrid$x
    squareY <- assocGrid$y
    squareData <- if (antennas) length(squares[[squareX]][[squareY]]$antennas)
      else length(squares[[squareX]][[squareY]]$people)
    if (antennas)
      squares[[squareX]][[squareY]]$antennas[[squareData + 1]] <- list(
        lon = inData$lon[row], lat = inData$lat[row],
        radio = inData$radio[row]
      )
    else
      squares[[squareX]][[squareY]]$people[[squareData + 1]] <- list(
        lon = inData$lon[row], lat = inData$lat[row]
      )
  }
  
  # Write info. to a JSON file
  gridObj$squares <- squares
  write(toJSON(gridObj, indent = 4, method = "C"), griddedJSON)
}


#' @description Match which people is assigned to each antenna inside a map
#' square.
#' @param peopleSquare list of people inside the square (lat, lon)
#' @param antennasSquare list of  antennas inside the square (lat, lon, radio)
#' @param regionFreqs list(UMTS, GSM, LTE) with the carrier freqs
#' @return the people-antennas matching as a
#' list(lat, lon, antennaLat, antennaLon, antennaRadio, distance, pathLoss)
#' @return NULL if thre is no person or antenna in the square
matchPeopleAndAntennas <- function(peopleSquare, antennasSquare, regionFreqs) {
  if (length(peopleSquare) < 1 | length(antennasSquare) < 1) {
    return(NULL)
  }
  
  for (p in 1:length(peopleSquare)) {
    minPathLoss <- Inf
    minAntenna <- Inf
    minDistance <- Inf
    
    for (a in 1:length(antennasSquare)) {
      fc <- get(antennasSquare[[a]]$radio, regionFreqs)
      d <- distance(lat1 = peopleSquare[[p]]$lat, lon1 = peopleSquare[[p]]$lon,
                lat2 = antennasSquare[[a]]$lat,
                lon2 = antennasSquare[[a]]$lon)$distance
      pathLoss <- COST231WalfischIkegamiLOS(fc = fc, distance = d)
      
      if (pathLoss < minPathLoss) {
        minPathLoss <- pathLoss
        minAntenna <- a
        minDistance <- d
      }
    }
    
    peopleSquare[[p]]$antennaRadio <- antennasSquare[[a]]$radio
    peopleSquare[[p]]$antennaLat <- antennasSquare[[a]]$lat
    peopleSquare[[p]]$antennaLon <- antennasSquare[[a]]$lon
    peopleSquare[[p]]$distance <- minDistance
    peopleSquare[[p]]$pathLoss <- minPathLoss
  }
  
  print("len(peopleSquare")
  print(length(peopleSquare))
  return(peopleSquare)
}

                                               
#' @description It assigns people inside grids to the antennas falling inside
#' that grid.
#' @param antennasGridJSON JSON file where antennas' locations are gridded
#' @param peopleGridJSON JSON file where people's locations are gridded
#' @param regionFreqs list(UMTS, GSM, LTE) with the carrier freqs
#' @note the peopleGridJSON is modified to include the associated antenna inside
#' each person's object.
#' @note if the people has already been assigned, then no assignment is
#' performed
assignAntennas <- function(antennasGridJSON, peopleGridJSON, regionFreqs) {
  griddedAntennas <- fromJSON(file = antennasGridJSON)
  griddedPeople <- fromJSON(file = peopleGridJSON)
  
  # Check if distances are already assigned
  if (griddedPeople$assignedDis) {
    print("The received peopeGridJSON has already assigned its distances")
    return()
  }
  
  # Check if the region and sampling longitudes match
  if (griddedAntennas$lonN != griddedPeople$lonN |
      griddedAntennas$latN != griddedPeople$latN) {
    stop("Non matching longitude and latitude gridding in people and atennas\n")
  }
  if (griddedAntennas$latB != griddedPeople$latB |
      griddedAntennas$latT != griddedPeople$latT |
      griddedAntennas$lonL != griddedPeople$lonL|
      griddedAntennas$lonR!= griddedPeople$lonR) {
    stop("Non matching latitude and longitude regions in people and antennas\n")
  }
  print("x")
  print(length(griddedAntennas$squares))
  print("y")
  print(length(griddedAntennas$squares[[1]]))
  for (x in 1:length(griddedAntennas$squares)) {
    for (y in 1:length(griddedAntennas$squares[[x]])) {
      antennasSquare <- griddedAntennas$squares[[x]][[y]]
      peopleSquare <- griddedPeople$squares[[x]][[y]]
      
      assignedPeopleSquare <- matchPeopleAndAntennas(peopleSquare$people,
                                                     antennasSquare$antennas,
                                                     regionFreqs)
      if (!is.null(assignedPeopleSquare))
        griddedPeople$squares[[x]][[y]]$people <- assignedPeopleSquare
    }
  }
  
  print("a imprimir")
  write(toJSON(griddedPeople, indent = 4, method = "C"),
        file = "/tmp/a.json")
  
  print("dump at:")
  print(paste(strsplit(peopleGridJSON, ".json")[[1]], "-assigned.json",
              sep = ""))
  write(toJSON(griddedPeople, indent = 4, method = "C"),
        paste(strsplit(peopleGridJSON, ".json")[[1]], "-assigned.json",
              sep = ""))
}


############## SMALL CELLS UTILITIES ##############

#' @description Calculates Prob(x <= markI(lon, lat)) whith x in
#' Ball((lon,lat),r)
#' @param lon ball center longitude coordinate
#' @param lat ball center latitude coordinate
#' @param r ball radius
#' @param intensityF intensity function used for markI
#' @param rhoStep step length as percentage of r for radius
#' @param thetaStep step length as percentage of 2*PI
#' @return a probability number
#' @note rhoStep and thetaStep must be specified as 1/n, where n is the number
#' of sampling steps to calculate marks inside the ball
cdfBallMarkI <- function(lon, lat, r, intensityF, rhoStep, thetaStep) {
  m <- markI(lon, lat, intensityF)
  
  sampleMarks <- c()
  for (rho_ in seq(from = 0, to = r, by = r*rhoStep)) {
    for (theta_ in seq(from = 0, to = 2*pi, by = 2*pi*thetaStep)) {
      sampleMarks <- c(sampleMarks, markI(lon + rho_*cos(theta_),
                                          lat + rho_*sin(theta_), intensityF))
    }
  }
  
  return(sum(sampleMarks <= m) / length(sampleMarks))
}


#' @description Approximates the intensity measure of a mattern II hard core
#' process with markI marks within the rectangle of specified limits.
#' @param lonL left longitude coordinate for the rectangle
#' @param lonR right longitude coordinate for the rectangle
#' @param latB bottom latitude coordinate for the rectangle
#' @param latT top coordinate for the rectangle
#' @param r inhibition radius
#' @param intensityF intensity function associated with the point process
#' @param rhoStep step length as percentage of r for radius (see cdfBallMarkI)
#' for more details
#' @param thetaStep step length as percentage of 2*PI
#' @return the Q value of the integral2 function
#' @note the not min probability is approximated as the product of the intensity
#' measure at a ball, multiplied by 1-F_m(u) where F_m is the marks CDF
matternIImapIintMapprox1 <- function(lonL, lonR, latB, latT, r, intensityF,
                                     rhoStep, thetaStep) {
  integrand2 <- function(lon, lat) {
    Fm <- cdfBallMarkI(lon, lat, r, intensityF, rhoStep, thetaStep)
    meanBall <- ballIntM(lon, lat, r, intensityF)
    
    probNotMin <- meanBall * Fm
    
    return((1 - probNotMin) * intensityF(lon, lat))
  }
  
  return(integral2(integrand2, xmin = lonL, xmax = lonR,
                          ymin = latB, ymax = latT)$Q)
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
#' @note the not min probability is approximated as the integral of the
#' intensity function along those points with lower marks inside a ball
matternIImapIintMapprox2 <- function(lonL, lonR, latB, latT, r, intensityF) {
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
    
    return((1 - probNotMin) * intensityF(lon, lat))
  }
  
  return(integral2(integrand4, xmin = lonL, xmax = lonR,
                          ymin = latB, ymax = latT)$Q)
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


#' @description Generates an intensity function f(lon, lat) for the small cell
#' antennas generation arround a given LTE antenna.
#' @param lteLon longitude coordinate of the LTE antenna
#' @param lteLat latitude coordinate of the LTE antenna
#' @param b b param of the centStepIntensity() 
#' @param c c param of the centStepIntensity() 
#' @param avgSmallCells average number of small cells to generate
#' @note this intensity function generator is based in the centStepIntensity()
genSmallCellManta <- function(lteLon, lteLat, b, c, avgSmallCells) {
  genManta <- local({
    centLon <- lteLon
    centLat <- lteLat
    a <- a
    b <- b
    c <- c
    avgPoints <- avgSmallCells
    intMeasure <- townIntM(centLon = lteLon, centLat = lteLat, b = b, c = c)
    
    function(lon, lat) {
      xs <- distanceMulti(lat1 = lat, lon1 = lon,
                             lat2 = centLat, lon2 = centLon)
      return(avgPoints / intMeasure * centStepIntensity(x = xs, b = b, c = c))
    }
  })
  
  return(genManta)
}


#' TODO - not finished its implementation
#' @description Generates an intensity function f(lon, lat) for the small cell
#' antennas generation arround a given LTE antenna. The generated function is
#' supposed to be used for a MaternII hard-core process with markI marks.
#' @param lteCenters dataframe with lon and lat coordinates of centers
#' @param lteLon longitude coordinate of the LTE antenna
#' @param lteLat latitude coordinate of the LTE antenna
#' @param b b param of the centStepIntensity()
#' @param c c param of the centStepIntensity()
#' @param r Matern II inhibition radius
#' @param avgSmallCells average number of small cells to generate 
#' @param approxM intensity measure approximation method:
#'        - "approx1" for matternIImapIintMapprox1
#'        - "approx2" for matternIImapIintMapprox2
#'        - "approx3" for matternIImapIintMapprox3
#' @param rhoStep step length as percentage of r for radius
#' @param thetaStep step length as percentage of 2*PI
#' @note this intensity function generator is based in the centStepIntensity()
#' @note if approxM is not correctly specified, approx3 is used
sgenMatIImapIsmallCellManta <- function(lteLon, lteLat, b, c, r, avgSmallCells,
                                       approxM, rhoStep = 1/10,
                                       thetaStep = 1/10) {
  # Antena intensity function
  lteIntF <- c(lteIntF, local({
    lteLon_ <- lteLon
    lteLat_ <- lteLat
    b_ <- b
    c_ <- c
    
    function(lon, lat) {
      dis <- distanceMulti(lat1 = lat, lon1 = lon,
                           lat2 = lteLat_, lon2 = lteLon_)
      return(centStepIntensity(x = dis, b = b_, c = c_))
    }
  }))
  
  # Obtain the intensity measure for the MatternII markI process
  lonL <- lteCenters[i,]$lon - (2*b + c)
  lonR <- lteCenters[i,]$lon + (2*b + c)
  latB <- lteCenters[i,]$lat - (2*b + c)
  latT <- lteCenters[i,]$lat + (2*b + c)
  if (approxM == "approx1") {
    lteIntMs <- c(lteIntMs, matternIImapIintMapprox1(lonL = lonL, lonR = lonR,
                                                     latB = latB, latT = latT,
                                                     r = r,
                                                     intensityF = lteIntF,
                                                     rhoStep = rhoStep,
                                                     thetaStep = thetaStep))
  } else if(approxM == "approx2") {
    lteIntMs <- c(lteIntMs, matternIImapIintMapprox2(lonL = lonL, lonR = lonR,
                                                     latB = latB, latT = latT,
                                                     r = r,
                                                     intensityF = lteIntF))
  } else {
    lteIntMs <- c(lteIntMs, matternIImapIintMapprox3(lonL = lonL, lonR = lonR,
                                                     latB = latB, latT = latT,
                                                     r = r,
                                                     intensityF = lteIntF))
  }
  
  genManta <- local({
    lteLon_ <- lteLon
    lteLat_ <- lteLat
    avgSmallCells_ <- avgSmallCells
    lteIntF_ <- lteIntF
    lteIntMs_ <- lteIntMs
    bs_ <- bs
    cs_ <- cs
    
    # Add the obtained intensity functions centered at each antenna
    function(lon, lat) {
      sumInts <- 0
      
      for(i in 1:nrow(lteCenters)) {
        lteCenter <- lteCenters_[i,]
        dis <- distanceMulti(lat1 = lat, lon1 = lon,
                            lat2 = lteCenter$lat, lon2 = lteCenter$lon)
        sumInts <- sumInts + avgSmallCells_[i] / lteIntMs_[i] *
          lteIntF_(x = dis, b = bs_[i], c = cs_[i])
      }
      
      return(sumInts)
    }
  })
  
  return(genManta)
}


#' @description Generates an intensity function f(lon, lat) for the small cell
#' antennas generation arround a given LTE antennas. The generated function is
#' supposed to be used for a MaternII hard-core process with markI marks.
#' @param lteCenters dataframe with lon and lat coordinates of centers
#' @param bs b param of the centStepIntensity() at each antenna
#' @param cs c param of the centStepIntensity() at each antenna
#' @param rs Matern II inhibition radius at each antenna
#' @param avgSmallCells average number of small cells to generate at each anten
#' @param approxM intensity measure approximation method:
#'        - "approx1" for matternIImapIintMapprox1
#'        - "approx2" for matternIImapIintMapprox2
#'        - "approx3" for matternIImapIintMapprox3
#' @param rhoStep step length as percentage of r for radius
#' @param thetaStep step length as percentage of 2*PI
#' @note this intensity function generator is based in the centStepIntensity()
#' @note if approxM is not correctly specified, approx3 is used
genMatIImapIsmallCellManta <- function(lteCenters, bs, cs, rs, avgSmallCells,
                                       approxM,
                                       rhoStep = 1/10, thetaStep = 1/10) {
  lteIntF <- c()
  lteIntMs <- c()
  for (i in 1:nrow(lteCenters)) {
    # Antena intensity function
    lteIntF <- c(lteIntF, local({
      lteCenter <- lteCenters[i,]
      b <- bs[i]
      c <- cs[i]
      
      function(lon, lat) {
        dis <- distanceMulti(lat1 = lat, lon1 = lon,
                             lat2 = lteCenter$lat, lon2 = lteCenter$lon)
        return(centStepIntensity(x = dis, b = b, c = c))
      }
    }))
    
    
    # Obtain the intensity measure for the MatternII markI process
    rect <- outerRect(centLon = lteCenters[i,]$lon,
                      centLat = lteCenters[i,]$lat, radius = 2*bs[i] + cs[i])
    lonL <- rect$lonL
    lonR <- rect$lonR
    latB <- rect$latB
    latT <- rect$latT
    if (approxM == "approx1") {
      lteIntMs <- c(lteIntMs, matternIImapIintMapprox1(lonL = lonL, lonR = lonR,
                                                       latB = latB, latT = latT,
                                                       r = rs[i],
                                                       intensityF =
                                                         lteIntF[[i]],
                                                       rhoStep = rhoStep,
                                                       thetaStep = thetaStep))
    } else if(approxM == "approx2") {
      lteIntMs <- c(lteIntMs, matternIImapIintMapprox2(lonL = lonL, lonR = lonR,
                                                       latB = latB, latT = latT,
                                                       r = rs[i],
                                                       intensityF =
                                                         lteIntF[[i]]))
    } else {
      lteIntMs <- c(lteIntMs, matternIImapIintMapprox3(lonL = lonL, lonR = lonR,
                                                       latB = latB, latT = latT,
                                                       r = rs[i],
                                                       intensityF =
                                                         lteIntF[[i]]))
    }
  }
  genManta <- local({
    lteCenters_ <- lteCenters
    avgSmallCells_ <- avgSmallCells
    lteIntF_ <- lteIntF
    lteIntMs_ <- lteIntMs
    bs_ <- bs
    cs_ <- cs
    
    # Add the obtained intensity functions centered at each antenna
    function(lon, lat) {
      sumInts <- 0
      
      for(i in 1:nrow(lteCenters)) {
        sumInts <- sumInts + avgSmallCells_[i] / lteIntMs_[i] *
          lteIntF_[[i]](lon, lat)
      }
      
      return(sumInts)
    }
  })
  
  return(genManta)
}


#' @description It reads a JSON file with all the regions' information, and
#' returns the information related to the regionName
#' @param regionsFile JSON file with the regions' information
#' @param regionName string identifying the region to be selected
#' @return selected region as a list
selectRegion <- function(regionsFile, regionName) {
  regions <- fromJSON(file = regionsFile)
  region <- NULL
  
  for (region_ in regions$regions) {
    if (region_$id == regionName) {
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
  
  return(region)
}


#' @description Finds the associated population object of the passed region
#' where the coordinates (lon, lat) are located
#' @param regionList region list object from selectRegion()
#' @param lon longitude coordinate
#' @param lat latitude coordinate
#' @return population list
findAssocPopulation <- function(regionList, lon, lat) {
  clostestIndex <- 1
  closestDistance <- Inf
  
  for (p in length(regionList$populations)) {
    coords <- getCoords(coordsStr = regionList$populations[[p]]$center)
    dis <- distance(lat1 = lat, lon1 = lon,
                    lat2 = coords$lat, lon2 = coords$lon)$distance
    
    if (dis < closestDistance) {
      closestDistance <- dis
      closestIndex <- p
    }
  }
  
  return(regionList$populations[[p]])
}


#' @description Appends the small cells generated by a point process to an
#' antenna data.frame
#' @param antennasDf data.frame with antennas' info from OpenCellID CSV
#' @param smallCellsPP spatstat point process object with the small cells'
#' coordinates as (x=lon, y=lat)
#' @param operatorNET the operator net code (by default is Telefonica's one)
#' @return the merged data.frame
appendSmallCells <- function(antennasDf, smallCellsPP, operatorNET = 7) {
  addSmallCell <- function(x){
    if(is.factor(x)) return(factor(x, levels=c(levels(x), "smallCell")))
    return(x)
  }
  
  if (smallCellsPP$n == 0) {
    return(antennasDf)
  }
  extendedDf <- antennasDf
  
  lonIndex <- which(names(antennasDf) == "lon")
  latIndex <- which(names(antennasDf) == "lat")
  radioIndex <- which(names(antennasDf) == "radio")
  netIndex <- which(names(antennasDf) == "net")
  
  # Check if the small cell level is present
  if (length(which(levels(extendedDf$radio) == "LTE")) == 0) {
    extendedDf <- as.data.frame(lapply(extendedDf, addSmallCell))
  }
  
  for (s in 1:smallCellsPP$n) {
    smallCellRow <- rep(x = "NA", times = ncol(antennasDf))
    smallCellRow[lonIndex] <- smallCellsPP$x[s]
    smallCellRow[latIndex] <- smallCellsPP$y[s]
    smallCellRow[radioIndex] <- "smallCell"
    smallCellRow[netIndex] <- operatorNET
    
    extendedDf <- rbind(extendedDf, smallCellRow)
  }
  
  return(extendedDf)
}


#' @description Line between coordinates a and b
#' @param x longitude coordinate
#' @param a list with $lon and $lat values
#' @param b list with $lon and $lat values
#' @return latitude value for x longitude to get a point (x,y) falling inside
#' the line joining a and b
coordsLine <- function(x, a, b) {
  return((x - a$lon) / (b$lon - a$lon) * (b$lat - a$lat) + a$lat)
}


#' @description Gives the road latitudes associated with the given points when
#' the road line is approximated as the lines joinning the list of points.
#' @param pointsLongitudes longitudes of the points
#' @param milestones breakpoints along the road as a list[[index]] of $lon $lat
#' @return a data.frame of $lon $lat as the points along the road
#' @return NULL if the pointLongitudes do not fall within the milestones
#' @note the milestones and pointsLongitudes parameters must be ordered by
#' longitude in ascending order
roadLatitudes <- function(pointsLongitudes, milestones) {
  nLons <- length(pointsLongitudes)
  nMilestones <- length(milestones)
  if (pointsLongitudes[1] < milestones[[1]]$lon |
      pointsLongitudes[nLons] > milestones[[nMilestones]]$lon) {
    return(NULL)
  }
  
  mappedLatitudes <- c()
  m <- 1 # milestone index
  for (i in 1:length(pointsLongitudes)) {
    pLon <- pointsLongitudes[i]
    
    # Find the next milestone if point is not between current ones
    while (pLon > milestones[[m + 1]]$lon & m <= length(milestones)) {
      m <- m + 1
    }
    
    mappedLatitudes <- c(mappedLatitudes,
                         coordsLine(x = pLon, a = milestones[[m]],
                                    b = milestones[[m + 1]]))
  }
  
  return(data.frame(lon = pointsLongitudes, lat = mappedLatitudes))
}

