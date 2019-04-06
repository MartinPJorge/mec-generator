library(pracma)
library(SDMTools)
library(rjson)



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


#' @description Rescales the values of a matrix to the interval (A,B)
#' @param toA A value of rescale interval
#' @param toB B value of rescale interval
#' @return the reescale matrix with values between A and B
rescaleRange <- function(matrix, toA, toB) {
  minM <- min(matrix)
  maxM <- max(matrix)
  
  return((matrix - minM) / (maxM - minM) * (toB - toA) + toA)
}


#' @description Rescales the values of a matrix to the interval (A,B)
#' @param toA A value of rescale interval
#' @param toB B value of rescale interval
#' @return the reescale matrix with values between A and B
rescaleRange <- function(matrix, toA, toB) {
  minM <- min(matrix)
  maxM <- max(matrix)
  
  return((matrix - minM) / (maxM - minM) * (toB - toA) + toA)
}

#' @param coordsStr "latitude,longitude"
#' @description Obtains the coordinates from string values stored in regions for JSONs
#' @return data.frame(lat, lon)
getCoords <- function(coordsStr) {
  coords <- strsplit(coordsStr, ",")[[1]]
  return(data.frame(lon = as.numeric(coords[2]), lat = as.numeric(coords[1])))
}


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


#' @description Centered version of the step intensity function stepIntensity()
#' @param x radius variable
#' @param b step width
#' @param c center section width
#' @return number
centStepIntensity <- function(x, b, c) {
  return(stepIntensity(x = x + b + c/2, a = 0, b = b, c = b + c, d = 2*b + c))
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


#' @description Finds the region within a JSON object
#' @param regions data frame with regions' details
#' @param regionId identifier of the region within the region_json
#' @return list with the found region
findRegion <- function(regions, regionId) {
  region <- NULL
  for (region_ in regions$regions) {
    if (region_$id == regionId) {
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


#' @description Creates the intensity matrix for the cell antennas generation
#' @param region list with the selected region properties
#' @param lonSamples number of longitude samples for the matrix
#' @param latSamples number of latitude samples for the matrix
#' @return list with the intensity matrix, the axis of longitude/latitude
#'         samples list(intenMatrix, lonAxis, latAxis)
createMatrix <- function(region, lonSamples, latSamples) {
  # Obtain the region population areas and their details
  townPopus <- c()
  townRads <- c()
  townDisps <- c()
  popuLons <- c()
  popuLats <- c()
  for (population in region$populations) {
    townPopus <- c(townPopus, population$population)
    townRads <- c(townRads, population$radius)
    coords <- getCoords(population$center)
    popuLons <- c(popuLons, coords$lon)
    popuLats <- c(popuLats, coords$lat)
    townDisps <- c(townDisps, population$dispersion)
  }
  townCents <- data.frame(lat = popuLats, lon = popuLons)
  
  
  # Generate the intensity matrix of people used for the lambda of people
  peopleF <- genPopuManta(townCents = townCents, townRads = townRads,
                          townPopus = townPopus, townDisps = townDisps)
  intMat <- genIntMatrix(latis = c(region$bl$lat, region$tr$lat),
                          longis = c(region$bl$lon, region$tr$lon),
                          latisLen = latSamples,
                          longisLen = lonSamples, intF = peopleF) 
  
  # Normalize the lambda matrix
  #rePopMat <- log(intMat$matrix + 1)
  rePopMat <- intMat$matrix
  rePopMat <- rescaleRange(matrix = rePopMat, toA = 1, toB = 1.07)
  
  # Aleatorize with small scale to prevent equal marks
  minRan <- min(rePopMat) / 1000000
  randMat <- replicate(ncol(rePopMat), runif(nrow(rePopMat), min = minRan,
                                             max = 2 * minRan))
  rePopMat <- rePopMat + randMat
  
  return(list(intenMatrix = rePopMat, lonAxis = intMat$lonAxis,
              latAxis = intMat$latAxis))
}


#' @description Generates the template intensity function matrix for the area.
#' That is, still a factor parameter must be multiplying for specific
#' generation purposes.
#' @param regions absolute path to the JSON file with the regions locations
#' and census, or list with already read JSON file
#' @param regionId identifier of the region within the region_json
#' @param samples number of longitude and latitude samples to grid the region
#' @return list with the intensity matrix, the axis of longitude/latitude
#'         samples list(intenMatrix, lonAxis, latAxis)
baseIntensity <- function(regions, regionId, samples) {
  regions_ <- regions
  if (typeof(regions) == "character") {
    regions_ <- fromJSON(file = regions)
  }
  region <- findRegion(regions = regions_, regionId = regionId)
  return(createMatrix(region = region, lonSamples = samples,
                      latSamples = samples))
}



