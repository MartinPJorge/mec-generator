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