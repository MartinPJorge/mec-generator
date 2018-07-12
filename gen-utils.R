library(pracma)

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


# Calculates the intensity measure of an inhomogeneous PPP at an observed region
# after it has been thinned with points lying inside a ball. The function
# obtains the intensity measure at the ball B((x,y), r), with (x,y) inside
# obs_win.
# x: ball center X coordinate
# y: ball center Y coordinate
# intensity: intensity function of the original PPP
# obs_win: owin object (only rectangles are supported)
#
# returns: the intensity measure of the ball B((x,y), r)
rad_th_intm <- function(x, y, r, intensity, obs_win) {
  
  # Intensity function after Matern I thinning
  rad_th_int <- function(x, y) {
    # Intensity measure of the ball B((x,y), r)
    ball_intm <- function(x, y, r, intensity) {
      # Function to be integrated in polar coordinates with origin (x0, y0)
      polar_shifted_intensity <- function(r, theta) {
        x0 <- x
        y0 <- y
        return(r * intensity(r*cos(theta) - x0, r*sin(theta) - y0))
      }
    
      return(integral2(polar_shifted_intensity, xmin = 0, xmax = r,
                ymin = 0, ymax = 2*pi))
    }
    
    mu_B <- ball_intm(x, y, r, intensity)$Q
    
    # return((1 - exp(-1 * mu_B)) * intensity(x, y))
    return(exp(-1 * mu_B) * intensity(x, y))
  } 
  
  return(integral2(rad_th_int,
            xmin = obs_win$xrange[1], xmax = obs_win$xrange[2],
            ymin = obs_win$yrange[1], ymax = obs_win$yrange[2]))
}

