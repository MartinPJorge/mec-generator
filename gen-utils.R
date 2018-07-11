# It generates an intensity function that is the sum of more intensity
# functions centered at certain locations.
#
# centers: double dimension array of centers
# centers_intensity: probability density function to place at each center
# ...: arguments for the probability density function
#
# return: intensity function that is the sum of the centered ones.
#         function(x, y)
gen_manta <- function(centers, centers_intensity, ...) {
  # Defines the intensity functions for each center
  center_ints <- list()
  for (i in 1:dim(centers)[1]) {
    center_ints[[i]] <- local({
      origin <- centers[i,]
      
      function(x, y) {
        x_shifted <- x - origin[1]
        y_shifted <- y - origin[2]
        
        r <- sqrt(x_shifted^2 + y_shifted^2)
        return(centers_intensity(r, ...))
      }
      
    })
  }
  
  # Creates the function that is sum of intensities
  sum_int <- function(x, y) {
    sum_prob <- 0
    for (center_int in center_ints) {
      sum_prob <- sum_prob + center_int(x, y)
    }
    return(sum_prob)
  }
  
  return(sum_int)
}


# Obtains the euclidean distance between two points
# p1: vector
# p2: vector
#
# return: number
eucl_dis <- function(p1, p2) {
  return(sqrt(sum((p1 - p2)^2)))
}


# Matern I thinning. If a point lies inside the ball B((x,y), r), then the
# probability of removal is 1, otherwise is 0.
# x: first coordinate
# y: second coordinate
# ...:
#    - r: radius of the ball
rad_thin <- function(x, y, ...) {
  points <- matrix(ncol = length(x), nrow = length(y))
  points[,1] <- x
  points[,2] <- y
  arguments <- list(...)
  r <- arguments[1][[1]]
  print(r)
  deletions <- rep(0, length(x))
  
  # Loop through every point
  for (i in seq(1, length(x))) {
    j <- 1
    found_in_rad <- FALSE
    
    # Checks if there is another point closer than 'r' to (x[i],y[i])
    while(!found_in_rad & j <= length(x)) {
      dist <- eucl_dis(c(x[i], y[i]), c(x[j], y[j]))
      if (dist != 0 & dist <= r) {
        found_in_rad <- TRUE
      }
      j <- j + 1
    }
    
    deletions[i] <- if(found_in_rad) 0 else 1
  }
  
  return(deletions)
}