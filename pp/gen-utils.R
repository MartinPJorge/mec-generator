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
# r: thinning radius
# intensity: intensity function of the original PPP
# obs_win: owin object (only rectangles are supported)
#
# returns: the intensity measure of the ball B((x,y), r)
rad_th_intm <- function(r, intensity, obs_win) {
  
  # Intensity function after radial thinning
  rad_th_int <- function(x, y) {
    # Intensity measure of the ball B((x,y), r)
    ball_intm <- function(x, y, r, intensity) {
      # Function to be integrated in polar coordinates with origin (x0, y0)
      polar_shifted_intensity <- function(r, theta) {
        x0 <- x
        y0 <- y
        return(r * intensity(r*cos(theta) + x0, r*sin(theta) + y0))
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


# Calculates the intensity measure of an inhomogeneous PPP at an observed ball
# after it has been thinned with points lying inside a ball. The function
# obtains the intensity measure at the ball B((x,y), R)
# r: thinning radius
# intensity: intensity function of the original PPP
# x: first coordinate of the center in the ball observed region
# y: second coordinate of the center in the ball observed region
# R: radius of the ball observed region
#
# returns: the intensity measure of the ball B((x,y), R)
rad_th_intmR <- function(r, intensity, x, y, R) {
  c_0 <- x
  c_1 <- y
  
  # Intensity function after radial thinning
  rad_th_intR <- function(rho, theta) {
    # Intensity measure of the ball B((x,y), r)
    ball_intm <- function(x, y, r, intensity) {
      # Function to be integrated in polar coordinates with origin (x0, y0)
      polar_shifted_intensity <- function(r, theta) {
        x0 <- x
        y0 <- y
        return(r * intensity(r*cos(theta) + x0, r*sin(theta) + y0))
      }
    
      return(integral2(polar_shifted_intensity, xmin = 0, xmax = r,
                ymin = 0, ymax = 2*pi))
    }
    
    mu_B <- ball_intm(c_0 + rho*cos(theta), c_1 + rho*sin(theta), r, intensity)$Q
    
    # return((1 - exp(-1 * mu_B)) * intensity(x, y))
    return(exp(-1 * mu_B) * intensity(c_0 + rho*cos(theta), c_1 + rho*sin(theta)))
  } 
  
  return(integral2(rad_th_intR, xmin=0, xmax=R, ymin=0, ymax=2*pi)$Q)
}




# Obtains the intensity measure on a ball B((x,y), r)
intensityMR <- function(x, y, r, intensity) {
  shifted_intensity <- function(rho, theta, ...) { # x, y, intensity
    argsL <- list(...)
    x <- argsL[1][[1]]
    y <- argsL[2][[1]]
    intensity <- argsL[3][[1]]
    
    return(intensity(x + rho*cos(theta), y + rho*sin(theta)) * rho)
  }
  
  return(integral2(shifted_intensity, xmin=0, xmax=r, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   x, y, intensity)$Q)
}


# Obtains the intensity measure on a disc((x,y), r, R)
intensityMDisc <- function(x, y, r, R, intensity) {
  shifted_intensity <- function(rho, theta, ...) { # x, y, intensity
    argsL <- list(...)
    x <- argsL[1][[1]]
    y <- argsL[2][[1]]
    intensity <- argsL[3][[1]]
    
    return(intensity(x + rho*cos(theta), y + rho*sin(theta)) * rho)
  }
  
  return(integral2(shifted_intensity, xmin=r, xmax=R, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   x, y, intensity)$Q)
}


nPointsProb <- function(x, y, r, intensity, numPoints) {
  nu <- intensityMR(x, y, r, intensity)
  return(nu^numPoints * exp(-1 * nu) / fact(numPoints))
}

# Obtains the probability of n points falling in the disc
# D((x,y), r, R)
nPointsProbDisc <- function(x, y, r, R, intensity, numPoints) {
  nu <- intensityMDisc(x, y, r, R, intensity)
  return(nu^numPoints * exp(-1 * nu) / fact(numPoints))
}


# Obtains the contact distribution at a point (x,y), i.e.:
#   F(x,y, r) = Prob( dist( (x,y), X ) <= r )
# It is the F function associated with the thinning processes expressed
# in rad_th_intmR and rad_th_intR when thinning radius is R
contact_rad_th <- function(r, intensity, x, y, R) {
  
  # Aux function to obtain the probability of a kill ocurrence at a point
  # (x + rho*cos(theta), y + rho*sin(theta))
  kill_prob <- function(rho, theta, ...) { # intensity, x, y, R
    #
    argsL <- list(...)
    intensity <- argsL[1][[1]]
    x <- argsL[2][[1]]
    y <- argsL[3][[1]]
    R <- argsL[4][[1]]
    
    x_0 <- x + rho*cos(theta)
    y_0 <- y + rho*sin(theta)
    zeroPprob <- nPointsProb(x_0, y_0, R, intensity, 0)
    onePprob <- nPointsProb(x_0, y_0, R, intensity, 1)
    sprintf("zeroPro: %.2f", zeroPprob)
    sprintf("onePro: %.2f", onePprob)
    return((1 - (zeroPprob+onePprob)) * rho)
  }
  
  return(integral2(kill_prob, xmin=0, xmax=r, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   intensity, x, y, R)$Q)
}



#### TRY with a contact distribution based on the paper approach
contact_radTh_paper <- function(x, y, rTh, rContact, intensity) {
  integrand <- function(rho, theta, ...) { # x, y rTh, intensity
    argsL <- list(...)
    x <- argsL[1][[1]]
    y <- argsL[2][[1]]
    rTh <- argsL[3][[1]]
    intensity <- argsL[4][[1]]
    
    x_0 <- x + rho*cos(theta)
    y_0 <- y + rho*sin(theta)
    
    prob0poins <- nPointsProb(x_0, y_0, rTh, intensity, 0)
    
    return(rho * intensity(x_0, y_0) * prob0poins)
  }
  
  return(1 - exp(-1 * integral2(integrand, xmin=0, xmax=r, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   x, y, rTh, intensity)$Q))
}

#### TRY with a contact distribution based on the paper approach
contact_rad_paper <- function(x, y, rContact, intensity) {
  integrand <- function(rho, theta, ...) { # x, y, intensity
    argsL <- list(...)
    x <- argsL[1][[1]]
    y <- argsL[2][[1]]
    intensity <- argsL[3][[1]]
    
    x_0 <- x + rho*cos(theta)
    y_0 <- y + rho*sin(theta)
    
    # prob0poins <- nPointsProb(x_0, y_0, rTh, intensity, 0)
    
    return(rho * intensity(x_0, y_0))
  }
  
  return(1 - exp(-1 * integral2(integrand, xmin=0, xmax=rContact, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   x, y, intensity)$Q))
}


# Probability that all point fallen in B(x,y,rContact) is not deleted in thin
non_deleted <- function(x, y, rContact, rTh, intensity) {
  shift_numPointsProb <- function(rho, theta, ...) { # r, intensity, numPoints
    argL <- list(...)
    r <- argL[1][[1]]
    intensity <- argL[2][[1]]
    numPoints <- argL[3][[1]]
    
    x_0 <- x + rho*cos(theta)
    y_0 <- y + rho*sin(theta)
    
    return(rho * nPointsProb(x_0,y_0, r, intensity, numPoints) )
  }
  
  return((1 - nPointsProb(x,y,rTh, intensity,0) - nPointsProb(x,y,rTh, intensity,1)))
  
  return(integral2(shift_numPointsProb, xmin=0, xmax=rContact, ymin=0, ymax=2*pi,
                   sector = FALSE,
                   reltol = 1e-6, abstol = 0, maxlist = 5000,
                   singular = FALSE, vectorized = TRUE,
                   rTh, intensity, 0)$Q)
}




# Lower bound aproximation of F(r) in the thinned process
low_f_aprox <- function(x, y, contactR, rTh, intensity, N) {
  low_aprox <- 0 # probability low bound of 1 point B(x,contactR)
  
  for (n in seq(N)) {
    oneInDisc <- nPointsProbDisc(x, y, contactR/N*(n-1), contactR/N*n,
                             intensity, 1)
    zeroInnerDisc <- 0
    if (n == 1) {
      zeroInnerDisc <- 1
    }
    else {
      zeroInnerDisc <- nPointsProbDisc(x, y, 0, contactR/N*(n-1),
                             intensity, 0)
    }
    
    outRad <- max(rTh, contactR - n/N*contactR)
    zeroInOut <- nPointsProbDisc(x, y, n/N*contactR, n/N*contactR + outRad,
                           intensity, 0)
    low_aprox <- low_aprox + zeroInnerDisc * oneInDisc * zeroInOut
  }
  
  return(low_aprox)
}



##### ANOTHER APROXIMATION TO P(B(x,r)=1) #####
f_ochio_low_aprox <- function(x, y, contactR, rTh, intensity, epsilon) {
  integrand <- function(rho, theta, ...) { # o_x, o_y, contactR, rTh, intensity,
                                           # epsilon
    argsL <- list(...)
    o_x <- argsL[1][[1]]
    o_y <- argsL[2][[1]]
    contactR <- argsL[3][[1]]
    rTh <- argsL[4][[1]]
    intensity <- argsL[5][[1]]
    epsilon <- argsL[6][[1]]
    
    x_0 <- o_x + rho*cos(theta)
    y_0 <- o_y + rho*sin(theta)
    orig_dist <- sqrt((o_x - x_0)^2 + (o_y - y_0)^2)
    
    # Obtain probability of ocurrence at (x,y)
    ocur_rad <- min(epsilon, contactR - orig_dist)
    pointProb <- nPointsProb(x_0, y_0, r=ocur_rad, intensity, 1)
    
    # No ocurrence arround point
    noArround <- nPointsProbDisc(x_0, y_0, ocur_rad, ocur_rad+r, intensity, 0)
    
    ###
    # No ocurrence at the rest of the ball
    outside_ball <- function(x_0, y_0, x, y, ocur_rad, rTh) {
      check <- 1
      
      # if(sqrt((x - x_0)^2 + (y - y_0)^2) <= (ocur_rad+rTh)) {
      #   check <- 0
      # }
      return(check)
    }
    
    intensity_out <- function(rho, theta, ...) { # x, y, contactR, rTh,
                                                 # intensity, epsilon, ocur_rad,
                                                 # o_x, o_y
      argsL <- list(...)
      x <- argsL[1][[1]] # mini-ball center
      y <- argsL[2][[1]]
      contactR <- argsL[3][[1]]
      rTh <- argsL[4][[1]]
      intensity <- argsL[5][[1]]
      epsilon <- argsL[6][[1]]
      ocur_rad <- argsL[7][[1]]
      o_x <- argsL[8][[1]]
      o_y <- argsL[9][[1]]
      
      x__0 <- o_x + rho*cos(theta)
      y__0 <- o_y + rho*sin(theta)
      
      is_outside_ball <- outside_ball(x, y, x__0, y__0, ocur_rad, rTh)
      return(intensity(x__0, y__0) * is_outside_ball * rho)
    }
    prob_no_outside <- integral2(intensity_out, xmin=0, xmax=contactR,
                                 ymin=0, ymax=2*pi,
                                 sector = FALSE,
                                 reltol = 1e-6, abstol = 0, maxlist = 5000,
                                 singular = FALSE, vectorized = TRUE,
                                 x_0, y_0, contactR, rTh, intensity, epsilon,
                                 ocur_rad, o_x, o_y)$Q
    
    #return(pointProb * noArround * prob_no_outside)
    return(pointProb * noArround)
  }
  
  
  return(integral2(integrand, xmin=0, xmax=contactR,
                                 ymin=0, ymax=2*pi,
                                 sector = FALSE,
                                 reltol = 1e-6, abstol = 0, maxlist = 5000,
                                 singular = FALSE, vectorized = TRUE,
                                 x, y, contactR, rTh, intensity, epsilon)$Q)
}
