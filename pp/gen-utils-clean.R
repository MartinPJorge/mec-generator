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
#' @return the average number of points in the rectangle
inhIntM <- function(lon1, lon2, lat1, lat2, intensityF) {
  return(integral2(intensityF, xmin = lon1, xmax = lon2,
                   ymin = lat1, ymax = lat2))
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
#' @param lat latitude coord
#' @param lon longitude coord
#' @param intF intensity function to evaluate (lat, lon) intensity
#' @return 1 / intF(lat, lon)
markI <- function(lat, lon, intF) {
  return(1 / intF(lat, lon))
}


#' @description MaternII thinning based on Vicenty's distance
#' @param X labeled point pattern object of spatstat
#' @return data.frame(lat, lon, marks, n)
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
        dist <- distanceMulti(look$lat, look$lon, X[j]$lat, X[j]$lon)
        if (dist < r) {
          foundMinor <- X[j]$marks < look$marks
        }
      }
      j <- j + 1
    }
    
    if (!foundMinor) {
      survivorsLat <- c(survivorsLat, look$lat)
      survivorsLon <- c(survivorsLon, look$lon)
      survivorsMarks <- c(survivorsMarks, look$marks)
    }
  }
  
  return(data.frame(x = survivorsLat, y = survivorsLon, marks = survivorsMarks,
                    n = length(survivorsLat)))
}


#' @description generates a MaternII process using the markI labeling.
#' Marks of points are higher for lower intensity. The process is for points in
#' maps, and it uses Vicenty's distance for the inhibition.
#' @param lambda intensity function
#' @param win owin rectangle with the limits
#' @param r MaternII inhibition radius
#' @return data.frame(x, y, marks, n)
jorgeMaternIImapI <- function(lambda, win, r) {
  P <- rpoispp(lambda, win = unitSquare, nsim = 1)
  
  # Generate the marks for the points
  marks <- c()
  for (i in 1:P$n) {
    marks <- c(marks, markI(P[i]$lat, P[i]$lon, lambda))
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


