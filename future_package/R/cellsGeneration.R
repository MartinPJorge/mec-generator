library(pracma)
library(SDMTools)
library(spatstat)


#' @description Generates cell antennas using given intensity matrix and squared
#' region factors
#' @param intenMatrix intensity matrix returned by baseIntensity()
#' @param lonAxis longitude axis returned by baseIntensity()  
#' @param latAxis latitude axis returned by baseIntensity()
#' @param squareFactors data.frame returned by intensityFactors(), with the
#' given square limits and intensity factors
#' @return data.frame with the longitude and latitude coordinates of the cell
#' antennas, plus additional information, i.e.,
#' data.frame(lon, lat, squareLonL, squareLonR, squareLatB, squareLatT,
#'            squareNum)
genCells <- function(intenMatrix, lonAxis, latAxis, squareFactors) {
  # Read the people lambda function and the longitude and latitude samples
  latT <- max(latAxis)
  latB <- min(latAxis)
  lonL <- min(lonAxis)
  lonR <- max(lonAxis)
  
  # Generate the interpolation function
  # print("Generating the people intensity function as input matrix interpolation")
  peopleIntpF <- genInterpInt(intenMatrix, lonAxis, latAxis)
  
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
    
    # cat(sprintf("Generating macro cells for square %d\n", row))
    # cat(sprintf("  lonL=%f lonR=%f latB=%f latT=%f\n", square$lonL, square$lonR,
    #         square$latB, square$latT))
    
    squareWin <- owin(xrange = c(square$lonL, square$lonR),
                      yrange = c(square$latB, square$latT))
    squareRepulsion <- 1 / (2 * square$sideSquares)
    
    lambdaSq <- function(lon, lat) {
      return(square$factor / square$avgCells * peopleIntpF(lon, lat))
    }
    
    # Generate the macro-cell antennas
    squareMacroCells <- jorgeMaternIImapI(lambda = lambdaSq, win = squareWin,
                                      r = squareRepulsion * 1000)
    
    #cat(sprintf("  %d macro-cells generated\n", squareMacroCells$n))
    
    # Append the antennas and square information
    macroLons <- c(macroLons, squareMacroCells$x)
    macroLats <- c(macroLats, squareMacroCells$y)
    lonLs <- c(lonLs, rep(square$lonL, squareMacroCells$n))
    lonRs <- c(lonRs, rep(square$lonR, squareMacroCells$n))
    latBs <- c(latBs, rep(square$latB, squareMacroCells$n))
    latTs <- c(latTs, rep(square$latT, squareMacroCells$n))
    squareNos <- c(squareNos, rep(row, squareMacroCells$n))
  }
    
    
  return(data.frame(lon = macroLons, lat = macroLats, squareLonL = lonLs,
                     squareLonR = lonRs, squareLatB = latBs, squareLatT = latTs,
                     squareNum = squareNos,
                     radio = rep("macro-cell", length(lonRs))))
}

