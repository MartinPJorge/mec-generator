library(pracma)
source("gen-utils-clean.R")
library(spatstat)

xs <- c(0, 10) 
ys <- xs
matLatSamples <- 50
matLonSamples <- matLatSamples
win <- owin(xrange = xs, yrange = ys)
factor <- 1
n <- 1
centers <- data.frame(lon = runif(n = n, min = 0, max = 10),
                      lat = runif(n = n, min = 0, max = 10))


###### TRY MULTIPLE density functions ####
# manta <- gen_manta(centers = centers, centers_intensity = dnorm,
#                       factor = factor, mean = 0, sd = 200000)
# manta <- gen_manta(centers = centers, centers_intensity = dunif,
#                       factor = factor, min = 0, max = 200000)
# manta <- gen_manta(centers = centers, centers_intensity = dgamma,
#                       factor = factor, shape = 1, scale = 200000)
manta <- gen_manta(centers = centers, factor = factor,
                        centers_intensity = stepIntensity,
                        a = 0, b = 1,
                        c = 200000, d = 200001)
#########################################


intMat <- genIntMatrix(latis = xs, longis = ys, latisLen = matLatSamples,
                       longisLen = matLonSamples, intF = manta)
intpIntF <- genInterpInt(intMat = intMat$matrix, lonAxis = intMat$lonAxis,
                         latAxis = intMat$latAxis)
rInhib <- 0.5

numPoints <- c()
for (i in 1:40) {
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


integrand4 <- function(lon, lat) {
  Fm <- intFCumBall(lon = lon, lat = lat, r = rInhib)
  meanBall <- ballIntM(lon, lat, r = rInhib)
  lambda <- intpIntF(lon, lat)
  
  filterBigger <- function(rho, theta) {
    pLambda <- intpIntF(lon + rho*cos(theta), lat + rho*sin(theta))
    pLambda <- replace(pLambda, pLambda >= lambda, 0)
    
    return(pLambda * rho)
  }
  
  probNotMin <- integral2(filterBigger, xmin = 0, xmax = rInhib,
                          ymin = 0, ymax = 2*pi)$Q
  
  return((1 - probNotMin) * lambda)
}

measureMat <- integral2(integrand4, xmin = xs[1], xmax = xs[2], ymin = ys[1],
                        ymax = ys[2])
print("Mean number of points")
print(mean(numPoints))
print("Approx")
print(measureMat$Q)



filled.contour(x = seq(xs[1], xs[2], length.out = nrow(intMat$matrix)),
               y = seq(ys[1], ys[2], length.out = ncol(intMat$matrix)),
               z = intMat$matrix)
points(x = points_$lon, y = points_$lat)
