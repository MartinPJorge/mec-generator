# TEST FILE

library(spatstat)
source("gen-utils-clean.R")

rThinning <- 2
centers <- data.frame(x = c(1, 5, 9), y = c(9, 1, 9))
gen_manta_ <- function(centers, centers_intensity, factor, ...) {
  # Defines the intensity functions for each center
  center_ints <- list()
  for (i in 1:dim(centers)[1]) {
    center_ints[[i]] <- local({
      origin <- centers[i,]
      
      function(x, y) {
        # Vicenty's distance
        r <- sqrt((x - origin$x)^2 + (y - origin$y)^2)
         
        return(centers_intensity(r, ...))
      }
      
    })
  }
  
  # Creates the function that is sum of intensities
  sum_int <- function(x, y) {
    sum_prob <- 0
    i <- 0
    for (center_int in center_ints) {
      sum_prob <- sum_prob + center_int(x, y)
      i <- i + 1
    }
    return(as.numeric(sum_prob * factor))
  }
  
  return(sum_int)
}
mantita <- gen_manta_(centers = centers, centers_intensity =  stepIntensity,
                     factor = 1, a = 0, b = 0.1, c = 1, d = 8)
intMatrix <- genIntMatrix(c(0,10), c(0,10), latisLen=50, longisLen=50,
                       intF = mantita)
  
unitSquare <- owin(xrange = c(0, 10), yrange = c(0, 10))


Pthin <- jorgeMaternIImapI(lambda = mantita, win = unitSquare, r = rThinning)
timestamp()
contour(x = intMatrix$lonAxis, y = intMatrix$latAxis,
        z = intMatrix$matrix)
graphics::points(x = Pthin$x, y = Pthin$y, col = "green")
