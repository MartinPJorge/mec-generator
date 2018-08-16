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
genInterpInt <- function(intMat, lonAxis, latAxis) {
  args <- local({
    list(intMat = intMat, lonAxis = lonAxis, latAxis = latAxis)
  })
  
  interpolator <- function(lon, lat) {
    # interp2 function states y=rows and x=cols, hence y=longitude x=latitude
    interpVal <- interp2(x = args$latAxis, y = args$lonAxis, Z = args$intMat,
            xp = lat, yp = lon, method = "linear")
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


#' @description Marks a point as the inverse of its intensity
#' @param lon longitude coord
#' @param lat latitude coord
#' @param intF intensity function to evaluate (lon, lat) intensity
#' @return 1 / intF(lon, lat)
markI <- function(lon, lat, intF) {
  mark <- 1 / intF(lon, lat)
  
  # Random number to diferentiate between coordinates with same intensity
  zeros <- abs(log(mark))
  rand_ <- runif(1, min = 1e-1*(zeros - 2), max = 1e-1*(zeros - 1))
  
  return(mark + rand_)
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
  print("jorgeMaternII - gener")
  timestamp()
  P <- rpoispp(lambda, win = win, nsim = 1)
  timestamp()
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
#' @radius circle radius in meters
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
