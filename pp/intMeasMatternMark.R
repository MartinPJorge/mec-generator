library(pracma)
source("gen-utils-clean.R")
library(spatstat)

xs <- c(0, 10) 
ys <- xs
win <- owin(xrange = xs, yrange = ys)
factor <- 100000
centers <- data.frame(lon = c(2.5, 7.5, 1, 2, 3), lat = c(2.5, 7.5, 1, 2, 3))
n <- 1
centers <- data.frame(lon = runif(n = n, min = 0, max = 10),
                      lat = runif(n = n, min = 0, max = 10))



manta <- gen_manta(centers = centers, centers_intensity = dnorm,
                      factor = factor, mean = 0, sd = 200000)
intMat <- genIntMatrix(latis = xs, longis = ys, latisLen = 50, longisLen = 50,
                       intF = manta)
intpIntF <- genInterpInt(intMat = intMat$matrix, lonAxis = intMat$lonAxis,
                         latAxis = intMat$latAxis)
rInhib <- 0.5

numPoints <- c()
for (i in 1:200) {
  points_ <- jorgeMaternIImapI(lambda = intpIntF, win = win, r = rInhib)
  numPoints <- c(numPoints, points_$n)
  print("itration")
  print(i)
}
print(mean(numPoints))

intFCum <- function(lon, lat) {
  m <- 1 / intpIntF(lon, lat)
  matValues <- as.vector(intMat$matrix)
  markValues <- 1 / matValues
  return(sum(markValues <= m) / length(markValues))
}


rhoStep <- 1/10
thetaStep <- 1/10
intFCumBall <- function(lon, lat, r) {
  m <- 1 / intpIntF(lon, lat)
  
  sampleMarks <- c()
  for (rho_ in seq(from = 0, to = r, by = r*rhoStep)) {
    for (theta_ in seq(from = 0, to = 2*pi, by = 2*pi*thetaStep)) {
      sampleMarks <- c(sampleMarks, 1 / intpIntF(lon + rho_*cos(theta_),
                                                 lat + rho_*sin(theta_)))
    }
  }
  
  return(sum(sampleMarks <= m) / length(sampleMarks))
}


# Calc intensity measure
ballIntM <- function(lon, lat, r) {
  shifed <- function(rho, theta) {
    return(rho * intpIntF(lon + rho*cos(theta), lat + rho*sin(theta)))
  }
  
  return(integral2(shifed, xmin = 0, xmax = r, ymin = 0, ymax = 2*pi)$Q)
}

integrand <- function(lon, lat) {
  Fm <- intpIntF(lon, lat)
  meanBall <- ballIntM(lon, lat, r = rInhib)
  return((1 - Fm) * intpIntF(lon, lat))
}

integrand2 <- function(lon, lat) {
  Fm <- intFCumBall(lon = lon, lat = lat, r = rInhib)
  meanBall <- ballIntM(lon, lat, r = rInhib)
  return((1 - meanBall * Fm) * intpIntF(lon, lat))
}

integrand3 <- function(lon, lat) {
  Fm <- intFCumBall(lon = lon, lat = lat, r = rInhib)
  meanBall <- ballIntM(lon, lat, r = rInhib)
  
  notMin <- 0
  for (n in 1:floor(meanBall)) {
    notMin <- notMin + Fm^n * (1 - Fm)^(meanBall - n)
  }
  
  return((1 - notMin) * intpIntF(lon, lat))
}

integrand4 <- function(lon, lat) {
  Fm <- intFCumBall(lon = lon, lat = lat, r = rInhib)
  meanBall <- ballIntM(lon, lat, r = rInhib)
  lambda <- intpIntF(lon, lat)
  
  filterBigger <- function(rho, theta) {
    pLambda <- intpIntF(lon + rho*cos(theta), lat + rho*sin(theta))
    if(pLambda >= lambda) {
      return(0)
    }
    return(pLambda * rho)
  }
  
  probNotMin <- integral2(filterBigger, xmin = 0, xmax = rInhib,
                          ymin = 0, ymax = 2*pi)$Q
  
  return((1 - probNotMin) * lambda)
}

measureMat <- integral2(integrand4, xmin = xs[1], xmax = xs[2], ymin = ys[1],
                        ymax = ys[2])
print(measureMat$Q)
print(points_$n)



filled.contour(x = seq(xs[1], xs[2], length.out = nrow(intMat$matrix)),
               y = seq(ys[1], ys[2], length.out = ncol(intMat$matrix)),
               z = intMat$matrix)
points(x = points_$lon, y = points_$lat)

