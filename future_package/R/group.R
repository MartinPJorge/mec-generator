####################################################################
####################### BASE FUNCTIONS #############################
####################################################################


#' @description Obtains the square centered at (lon, lat) inside a region.
#' @param lon longitude coordinate
#' @param lat latitude coordinate
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param halfSide how long is in meters half the side of the square
#' @return list(left, right, bottom, up) with the limiting coordinates
#' @note it uses Vicenty's distance
aroundSquare <- function(lon, lat, lonL, lonR, latB, latT, halfSide) {
  # Get the sides' limits arround the MEC location (maxInt)
  leftLim <- destination(lat = lat, lon = lon,
                         distance = halfSide, bearing = 270)$lon2
  rightLim <- destination(lat = lat, lon = lon,
                          distance = halfSide, bearing = 90)$lon2
  topLim <- destination(lat = lat, lon = lon,
                        distance = halfSide, bearing = 0)$lat2
  bottomLim <- destination(lat = lat, lon = lon,
                           distance = halfSide, bearing = 180)$lat2
  
  # Check the limits are inside the region
  leftLim <- max(leftLim, lonL)
  rightLim <- min(rightLim, lonR)
  bottomLim <- max(bottomLim, latB)
  topLim <- min(topLim, latT)

  return(list(left = leftLim, right = rightLim,
              bottom = bottomLim, top = topLim))
}


#' @description Obtains the square centered at (lon, lat) inside a region. It
#' uses an estimation of the ammount of latitude and longitude degrees to
#' achieve a meter.
#' @param lon longitude coordinate
#' @param lat latitude coordinate
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param halfSide how long is in meters half the side of the square
#' @param lonMet longitude degrees per meter (deg/meter)
#' @param latMet latitude degrees per meter (deg/meter)
#' @return list(left, right, bottom, up) with the limiting coordinates
aroundSquareEstim <- function(lon, lat, lonL, lonR, latB, latT, halfSide,
                              lonMet, latMet) {
  
  # Square limits
  leftLim <- lon - abs(halfSide * lonMet)
  rightLim <- lon + abs(halfSide * lonMet)
  bottomLim <- lat - abs(halfSide * latMet)
  topLim <- lat + abs(halfSide * latMet)
  
  # Check the limits are inside the region
  leftLim <- max(leftLim, lonL)
  rightLim <- min(rightLim, lonR)
  bottomLim <- max(bottomLim, latB)
  topLim <- min(topLim, latT)

  return(list(left = leftLim, right = rightLim,
              bottom = bottomLim, top = topLim))
}


#' @description It tells what antennas fall inside the square.
#' @param antLons longitudes of the antennas
#' @param antLats latitudes of the antennas
#' @param checkAntenna vector of booleans telling if antenna must be considered
#' @param cLon longitude coordinate of the square center
#' @param cLat latitude coordinate of the square center
#' @param halfSide half side of the square in meters
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
antennasInSquare <- function(antLons, antLats, checkAntenna,
                             cLon, cLat, halfSide, lonL, lonR, latB, latT) {
  # Get the square surounding (cLon, cLat)
  arqSq <- aroundSquare(lon = cLon, lat = cLat,
                       lonL = lonL, lonR = lonR, latB = latB, latT = latT,
                       halfSide = halfSide)
  leftLim <- arqSq$left
  rightLim <- arqSq$right
  bottomLim <- arqSq$bottom
  topLim <- arqSq$top
  
  # Get the indexes of antennas inside the MEC square which are not covered
  # yet
  antInd <- (antLons > leftLim & antLons < rightLim) &
    (antLats > bottomLim & antLats < topLim) &
    checkAntenna
  antInd <- which(TRUE == antInd)
  
  return(antInd)
}


#' @description It obtains the closest antennas that can reach (cLon, cLat)
#' @param antLons longitudes of the antennas
#' @param antLats latitudes of the antennas
#' @param antRadios radio technologies ("LTE", "macro-cell", "femto-cell")
#' @param antDiss list with maximum distance that the antenna can be away from
#' (cLon, cLat) -> list(LTE, macroCell, femtoCell)
#' @param checkAntenna vector of booleans telling if antenna must be considered
#' @param cLon longitude coordinate of the square center
#' @param cLat latitude coordinate of the square center
#' @param r circle radius in meters. Inside that circles antennas are checked.
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param lonMeter longitudes degrees corresponding to a meter
#' @param latMeter latitude degrees corresponding to a meter
#' @param N number of closest antennas
#' @return antenna indexes that can reach (cLon, cLat)
#' @note the function uses Manhattan distance
#' @note if N is not specified, then the function returns all the antennas
#' indexes ordered increasingly by their distance to (cLon, cLat)
reachableAntennas <- function(antLons, antLats, antRadios, antDiss,
                              checkAntenna, cLon, cLat, r,
                              lonL, lonR, latB, latT,
                              lonMeter, latMeter, N = NULL) {
  
  # Indexes of antennas inside square
  antInds <- antennasInSquare(antLons, antLats, checkAntenna,
                              cLon, cLat, halfSide = r, lonL, lonR, latB, latT)
  cat(sprintf("     there are %d antennas inside square\n", length(antInds)))
  reachAntInds <- c()
  dis2Center <- c()
  
  for (ai in antInds) {
    ant2MecDis <- abs(antLons[ai] - cLon) / lonMeter +
      abs(antLats[ai] - cLat) / latMeter # Manhattan distance
    
    # Select antenna's maximum distance
    mDis <- NULL
    if (antRadios[ai] == "LTE") {
      mDis <- antDiss$LTE
    } else if(antRadios[ai] == "femto-cell") {
      mDis <- antDiss$femtoCell
    } else if(antRadios[ai] == "macro-cell") {
      mDis <- antDiss$macroCell
    }
    
    # Adds the antenna index if it reaches (cLon, cLat)
    if (ant2MecDis <= mDis) {
      reachAntInds <- c(reachAntInds, ai)
      dis2Center <- c(dis2Center, ant2MecDis)
    }
  }
  cat(sprintf("     there are %d reachable antennas inside square\n",
              length(reachAntInds)))
  
  # Order the indexes by its distance, and get the N closest ones
  antDf <- data.frame(inds = reachAntInds, dis = dis2Center)
  antDf <- antDf[order(antDf$dis),]
  if (is.null(N)) {
    N <- nrow(antDf)
  } else if (N > nrow(antDf)) {
    N <- nrow(antDf)
  }
  antDf <- antDf[1:N,]
  cat(sprintf("     only the closest %d are selected\n", N))
  
  return(antDf$inds)
}


#' @description It performs operations in the matrix entries corresponding to
#' locations falling inside the ball of center (cLon, cLat) and radius 'r'
#' @param matrix matrix indexed by [longitude, latitude]
#' @param operation string indicating the operation to be performed:
#'        zero, plus-one, minus-one
#' @param cLon longitude coordinate of the square center
#' @param cLat latitude coordinate of the square center
#' @param r radius of the ball centered at (cLon, cLat) where matrix operations
#' are performed. Units are meters
#' @param lonAxis vector with the longitude coordinates of the matrix
#' @param latAxis vector with the latitude coordinates of the matrix
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param lonMeter amount of longitude degrees corresponding to one meter
#' @param latMeter amount of latitude degrees corresponding to one meter
operateAroundCircle <- function(matrix, operation, cLon, cLat, r,
                                lonAxis, latAxis, lonL, lonR, latB, latT,
                                lonMeter, latMeter) {
  # Get the square surounding the antenna
  arqSq <- aroundSquare(lon = cLon, lat = cLat,
                       lonL = lonL, lonR = lonR, latB = latB, latT = latT,
                       halfSide = r)
  leftLim <- arqSq$left
  rightLim <- arqSq$right
  bottomLim <- arqSq$bottom
  topLim <- arqSq$top
  
  # Get the limiting square indexes inside the matrix
  lonInds <- which(lonAxis > leftLim & lonAxis < rightLim)
  latInds <- which(latAxis > bottomLim & latAxis < topLim)
  
  # Iterate through coords inside limiting square to get the ones inside
  # (cLon, cLat) circle
  for (lonInd in lonInds) {
    for (latInd in latInds) {
      manhattanDis <- abs(cLon - lonAxis[lonInd]) / lonMeter +
        abs(cLat - latAxis[latInd]) / latMeter
      
      # Add a unit to the coordinate
      if (manhattanDis <= r) {
        if (operation == "plus-one") {
          matrix[lonInd, latInd] <- matrix[lonInd, latInd] + 1
        }
        else if (operation == "minus-one") {
          matrix[lonInd, latInd] <- matrix[lonInd, latInd] - 1
        }
        else if (operation == "zero") {
          matrix[lonInd, latInd] <- 0
        }
      }
    }
  }
    
  return(matrix)
}


#' @description Given some antennas' location, it obtains the intensity function
#' of where MEC PoPs can be located using Manhattan distance.
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param lonSamples number of samples for the longitude axis in the matrix
#' @param latSamples number of samples for the latitude axis in the matrix
#' @param antLons longitudes of the antennas
#' @param antLats latitudes of the antennas
#' @param antRadios radio technology of each antenna (LTE, macro-cell, or
#' femto-cell)
#' @param maxDiss maximum distance where a MEC server can be away from each
#' antenna depending on its technology: list(LTE=10, macroCell=10, femtoCell=1)
#' @return list(matrix, latAxis, lonAxis)
#' @note the function estimates on the top left of the region how many latitude
#' and coordinate units are needed for a meter. This estimation is used all over
#' the map. The reasoning is to boost up performace avoiding Vicenty's distance.
#' Hence this function does not generate an exact intensity matrix
mecIntManhattan <- function(lonL, lonR, latB, latT, lonSamples, latSamples,
                            antLons, antLats, antRadios, maxDiss) {
  # Estimation of longitude and latitude distance of one meter
  lonD <- destination(lat = latT, lon = lonL, distance = 1, bearing = 90)
  lonMeter <- abs(lonD$lon2 - lonD$lon1)
  latD <- destination(lat = latT, lon = lonL, distance = 1, bearing = 180)
  latMeter <- abs(latD$lat2 - latD$lat1)
  
  # Create the intensity 
  lonAxis <- seq(from = lonL, to = lonR, length.out = lonSamples)
  latAxis <- seq(from = latB, to = latT, length.out = latSamples)
  mecMatrix <- matrix(data = 0, nrow = lonSamples, ncol = latSamples)
  D_ANT_INCREASES <- rep(x = 0, length(antLats))
  
  for (i in 1:length(antLons)) {
    print(i)
    # Select antenna's maximum distance
    mDis <- NULL
    if (antRadios[i] == "LTE") {
      mDis <- maxDiss$LTE
    } else if(antRadios[i] == "femto-cell") {
      mDis <- maxDiss$femtoCell
    } else if(antRadios[i] == "macro-cell") {
      mDis <- maxDiss$macroCell
    }
    
    mecMatrix <- operateAroundCircle(matrix = mecMatrix, operation = "plus-one",
                        cLon = antLons[i], cLat = antLats[i], r = mDis,
                        lonAxis = lonAxis, latAxis = latAxis,
                        lonL = lonL, lonR = lonR, latB = latB, latT = latT,
                        lonMeter = lonMeter, latMeter = latMeter)
  }
    
  
  return(list(matrix = mecMatrix, latAxis = latAxis, lonAxis = lonAxis,
              D_INCREASES = D_ANT_INCREASES))
}


#' @description It determines the location of MEC PoPs by not allowing them to
#' be closer than min(maxDiss) by digging down to 0 the mecIntMatrix around the
#' MEC PoP
#' @param mecIntMatrix intensity matrix [lon, lat] to generate the MEC PoP at a
#' certain coordinate
#' @param lonAxis vector with the longitude coordinates of the matrix
#' @param latAxis vector with the latitude coordinates of the matrix
#' @param maxDiss maximum distance where a MEC server can be away from each
#' antenna depending on its technology: list(LTE=10, macroCell=10, femtoCell=1)
#' @param antLons longitudes of the antennas
#' @param antLats latitudes of the antennas
#' @param antRadios radio technology of each antenna (LTE, macro-cell, or
#' femto-cell)
#' @param lonL left longitude of the region
#' @param lonR right longitude of the region
#' @param latB bottom latitude of the region
#' @param latT top latitude of the region
#' @param letNoAssign (optional) boolean determining if it is allowed that an
#' AAU is not assigned. If not specified, antennas can be unusigned.
#' @param numMECs (optional) number of MEC PoPs to generate. If not specified,
#' the function generates MEC locations until all antennas are covered.
#' femto-cell)
#' @return list(pos=data.frame(lon, lat, intensityMeas, coveredAs),
#' antennas=data.frame(lon, lat, MEC)
#' @note ret$pos$intMeas is the integral arround the MEC location
#'       ret$antennas$MEC references the MEC's row in ret$pos
#' @note some of the antennas can be uncovered, these ones will have
#'       return$antennas$MEC=-1. The parameter letNoAssign controls this.
#' @note Algorithm overview:
#' 1. while mecPoPs < numMECs:
#' 2.   nextMEC <- max_(x,y) mecIntMatix(x,y)
#' 3.   for antenna in Ball(nextMEC, max(maxDiss)):
#' 4.     mecIntMatrix[Ball(antenna,maxDiss(antenna))] -= 1
#' 5.   mecIntMatrix[Ball(nextMEC,min(maxDiss))] = 0
#' 6.   if max(mecIntMat) = 0:
#' 7.     for each nonCovered in antennas:
#' 8.       mecIntMatrix[Ball(nonCovered,maxDiss(nonCovered)))] += 1
mecLocationDig <- function(mecIntMatrix, lonAxis, latAxis, maxDiss,
                           antLons, antLats, antRadios,
                           lonL, lonR, latB, latT, letNoAssign = TRUE,
                           numMECs = NULL) {
  if (!is.null(numMECs) & letNoAssign == FALSE) {
    cat("You are invoking mecLocationDig() with numMECS specified and",
        "letNoAssign = FALSE. The function can't guarantee that, EXIT!")
    return(NULL)
  }
  
  # Estimation of longitude and latitude distance of one meter
  lonD <- destination(lat = latT, lon = lonL, distance = 1, bearing = 90)
  lonMeter <- abs(lonD$lon2 - lonD$lon1)
  latD <- destination(lat = latT, lon = lonL, distance = 1, bearing = 180)
  latMeter <- abs(latD$lat2 - latD$lat1)
  
  mecLons <- c()
  mecLats <- c()
  mecInts <- c()
  mecCoveredAs <- c() # num antennas each MEC covers
  antMecAsoc <- rep(x = -1, length(antLons))
  antCovered <- rep(x = FALSE, length(antLons))
  maxAntRad <- max(unlist(maxDiss))
  minAntRad <- min(unlist(maxDiss))
  putMoreMECs <- TRUE
  moreRounds <- ifelse(test = letNoAssign, yes = FALSE, no = TRUE)
  
  
  
  while (putMoreMECs) {
    maxMecInts <- max(mecIntMatrix)
    mecCoord <- which(maxMecInts == mecIntMatrix, arr.ind = TRUE)
    mecCoord_ <- which(maxMecInts == mecIntMatrix, arr.ind = TRUE)
    mecCoord <- c(lonAxis[mecCoord[1,1]], latAxis[mecCoord[1,2]])
    mecLons <- c(mecLons, mecCoord[1])
    mecLats <- c(mecLats, mecCoord[2])
    mecCoveredAs <- c(mecCoveredAs, 0)
    mecInts <- c(mecInts, maxMecInts)
    cat(sprintf("MEC number: %d\n", length(mecLons)))
    cat(sprintf("  max=%d, amount=%d\n", maxMecInts, nrow(mecCoord_)))
    cat(sprintf("  uncovered antennas: %d\n", sum(!antCovered)))
    cat(sprintf("  location: (%f,%f)\n", mecCoord[1], mecCoord[2]))
    cat(sprintf("  location indexes: (%d,%d)\n", mecCoord_[1], mecCoord_[2]))
    
    # Get the indexes of antennas inside square surounding the MEC
    antInd <- antennasInSquare(antLons = antLons, antLats = antLats,
                     checkAntenna = !antCovered,
                     cLon = mecCoord[1], cLat = mecCoord[2],
                     halfSide = maxAntRad, lonL = lonL, lonR = lonR,
                     latB = latB, latT = latT)
    
    # Iterate through the antennas covered by the MEC
    for (ai in antInd) {
      ant2MecDis <- abs(antLons[ai] - mecCoord[1]) / lonMeter +
        abs(antLats[ai] - mecCoord[2]) / latMeter # Manhattan distance
      
      # Select antenna's maximum distance
      mDis <- NULL
      if (antRadios[ai] == "LTE") {
        mDis <- maxDiss$LTE
      } else if(antRadios[ai] == "femto-cell") {
        mDis <- maxDiss$femtoCell
      } else if(antRadios[ai] == "macro-cell") {
        mDis <- maxDiss$macroCell
      }
      
      # If the antenna can access the MEC, associate and decrease its
      # suroundings
      if (ant2MecDis <= mDis) {
        antCovered[ai] <- TRUE
        antMecAsoc[ai] <- length(mecLons)
        mecCoveredAs[length(mecLons)] <- mecCoveredAs[length(mecLons)] + 1
        
        mecIntMatrix <- operateAroundCircle(matrix = mecIntMatrix,
                          operation = "minus-one", cLon = antLons[ai],
                          cLat = antLats[ai], r = mDis, lonAxis = lonAxis,
                          latAxis = latAxis, lonL = lonL, lonR = lonR,
                          latB = latB, latT = latT, lonMeter = lonMeter,
                          latMeter = latMeter)
      }
      
    }
    
    cat(sprintf("  num covered antennas: %d\n", mecCoveredAs[length(mecLons)] ))
    
    # Put to zero the located MEC PoPs suroundings
    mecIntMatrix <- operateAroundCircle(matrix = mecIntMatrix,
                      operation = "zero", cLon = mecCoord[1],
                      cLat = mecCoord[2], r = minAntRad, lonAxis = lonAxis,
                      latAxis = latAxis, lonL = lonL, lonR = lonR,
                      latB = latB, latT = latT, lonMeter = lonMeter,
                      latMeter = latMeter)
    mecIntMatrix[mecIntMatrix <= 0] <- 0
    
    
    # Decide if more MECs must be put
    if (!is.null(numMECs)) {
      putMoreMECs <- length(mecLons) < numMECs
    } else {
      if (letNoAssign) {
        putMoreMECs <- maxMecInts > 0
      } else {
        putMoreMECs <- sum(!antCovered) > 0
        
        # If all AAUs must be covered but the matrix max is zero, add 1 by
        if (maxMecInts < 1) {
          cat(
            sprintf("\n==> ANOTHER ITERATION TO COVER NON-COVERED AAUs <==\n"))
          for (ai in which(antCovered == FALSE)) {
            cat(sprintf("  add 1 around AAU(%f,%f)\n", antLons[ai],
                        antLats[ai]))
            # Select antenna's maximum distance
            mDis <- NULL
            if (antRadios[ai] == "LTE") {
              mDis <- maxDiss$LTE
            } else if(antRadios[ai] == "femto-cell") {
              mDis <- maxDiss$femtoCell
            } else if(antRadios[ai] == "macro-cell") {
              mDis <- maxDiss$macroCell
            }
            
            # Increase by one the antena surroundings
            mecIntMatrix <- operateAroundCircle(matrix = mecIntMatrix,
                              operation = "plus-one", cLon = antLons[ai],
                              cLat = antLats[ai], r = mDis, lonAxis = lonAxis,
                              latAxis = latAxis, lonL = lonL, lonR = lonR,
                              latB = latB, latT = latT, lonMeter = lonMeter,
                              latMeter = latMeter)
          } 
          
        }
      }
    }
      
  }
  
  cat("Uncovered antennas:")
  for (ai in which(antCovered == FALSE)) {
    cat(sprintf("  A(%f,%f)\n",
                antLons[ai], antLats[ai]))
  }
  
  
  # Pack the MEC locations and antenna associations
  mecs <- data.frame(lon = mecLons, lat = mecLats, intensityMeas = mecInts,
                     coveredAs = mecCoveredAs)
  antennasAssoc <- data.frame(lon = antLons, lat = antLats, MEC = antMecAsoc)
  
  return(list(pos = mecs, antennas = antennasAssoc, modMat = mecIntMatrix))
}



#############################################################
############# START FROM HERE ON THE ADAPTATION #############
#############################################################




#' @description It groups all elements with coordinates (lats, lons) within
#' groups of a maximum of groupN elements that are <maxDis from the group
#' center. It uses the Manhattan distance
groupElems <- function(lats, lons, lonAxis, latAxis, groupN, maxDis) {
  
}










