rm(list = ls())
source("gen-utils.R")
library(spatstat)
library(graphics) 
  
  ### GENERATION OF POINTS ###
graphics.off()
antenas <- matrix(ncol = 2, nrow = 2)
antenas[1,] <- c(-5, 5)
antenas[2,] <- c(5, 5)
mantaca <- gen_manta(antenas, dgamma, shape=2, scale=2)
xtics <- seq(-10, 10, length=100)
ytics <- seq(-10, 10, length=100)
values <- matrix(ncol=100, nrow=100)
i = 1;
for (x in xtics) {
  j = 1;
  for (y in ytics) {
    values[i, j] <- mantaca(x, y)
    j = j + 1
  }
  i = i + 1
}

r = 2 # thinning radius
rectangulo <- owin(xrange=c(-10,10), yrange=c(-10,10))
puntos <- rpoispp(mantaca, win = rectangulo)
rect_im <- im(values, xcol=xtics, yrow=ytics)
#contour.im(rect_im, add=FALSE, axes=TRUE, clipwin = rectangulo)
#contour(rect_im, add=TRUE, axes=TRUE, clipwin = rectangulo)
puntos_th <- rthin(X=puntos, P=rad_thin, r=r, nsim=1, drop=TRUE)
#puntos_th <- rMaternI(kappa = mantaca, r=r, win = rectangulo)
filled.contour(values, x = xtics, y = ytics, col=gray.colors(41),
               plot.axes = { points(puntos); axis(1); axis(2) }
)
filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
               plot.axes = { points(puntos_th); axis(1); axis(2) }
)
# points(antenas, pch=17)
# points(puntos, pch=19)
rectangulo_red <- owin(xrange=c(-10+r,10-r), yrange=c(-10+r,10-r))
int_th <- rad_th_intm(r=r, intensity=mantaca, obs_win=rectangulo)$Q
int_orig <- integral2(mantaca,
                 xmin = rectangulo$xrange[1], xmax = rectangulo$xrange[2],
                 ymin = rectangulo$yrange[1], ymax = rectangulo$yrange[2])$Q


cat("mu(R) = ", int_orig, "\n")
cat("mu_th(R) = ", int_th, "\n")
cat("mu(R) - mu_th(R) = ", int_orig - int_th, "\n")
cat("Realization X in R has ", length(puntos$x), " points\n")
cat("Realization X_th in R has ", length(puntos_th$x), " points\n")



#####
# Check if contact function matches
#####
numSims=100
fallIn = 0
numFallIn = array(0, dim=numSims)
x_0 = 0
y_0 = 0
contactR <- 2
for (sim in seq(numSims)) {
  puntos <- rpoispp(mantaca, win = rectangulo)
  filled.contour(values, x = xtics, y = ytics, col=gray.colors(41),
                 plot.axes = { points(puntos); axis(1); axis(2) }
  )
  puntos_th <- rthin(X=puntos, P=rad_thin, r=r, nsim=1, drop=TRUE)
  filled.contour(values, x = xtics, y = ytics, col=gray.colors(41),
                 plot.axes = { points(puntos_th); axis(1); axis(2) }
  )
  # puntos_th <- puntos
  
  if (length(puntos_th$x) != 0) {
    # For the ocurrence of falling inside radius of contact
    for (i in seq(length(puntos_th$x))) {
      x <- puntos_th$x[i]
      y <- puntos_th$y[i]
      
      distance <- sqrt((x-x_0)^2 + (y-y_0)^2)
      if (distance <= contactR) {
        fallIn = fallIn + 1
        break
      }
    }
    
    # To get avg number of points
    for (i in seq(length(puntos_th$x))) {
      x <- puntos_th$x[i]
      y <- puntos_th$y[i]
      distance <- sqrt((x-x_0)^2 + (y-y_0)^2)
      if (distance <= contactR) {
        numFallIn[sim] = numFallIn[sim] + 1
      }
    }
    
  }
}

# Intensity measure original process
cerca_antena <- intensityMR(x=x_0, y=y_0, r=contactR, intensity=mantaca)
avgCerca <- mean(numFallIn)
sprintf("Expected original at B(x=%d,y=%d, R=%.2f) = %.2f", x_0, y_0, contactR, cerca_antena)
sprintf("  experimental: %.2f", avgCerca)



contactProb <- contact_rad_th(r=contactR, intensity=mantaca, x=x_0, y=y_0, R=r)
contactProbPaper <- contact_radTh_paper(x=x_0, y=y_0, rTh=r, rContact=contactR,
                                        intensity = mantaca)
contactProbPaperOri <- contact_rad_paper(x=x_0, y=y_0, rContact=contactR,
                                        intensity = mantaca)
contactLowAprox <- low_f_aprox(x=x_0, y=y_0, contactR=contactR, rTh=r,
                               intensity=mantaca, N=10)
ochioLowAprox <- f_ochio_low_aprox(x=x_0, y=y_0, contactR=contactR, rTh=r,
                               intensity=mantaca, epsilon=1)

f_aprox <- 1 - (nPointsProb(x=x_0, y=y_0, r=contactR, intensity=mantaca, numPoints=0) +  0.7)


fallenRate <- fallIn / numSims


sprintf("F(x=%d,y=%d,r=%d)::\n", x_0, y_0, r)
sprintf("  empirical rate in %d sims -> F=%.2f", numSims, fallenRate)
sprintf("  paper approach / 2-> F=%.2f", contactProbPaper / 2)
sprintf("  lower bound -> F_=%.10f", contactLowAprox)
sprintf("  ochio lower bound -> F_=%.10f", ochioLowAprox)

aaa <- rad_th_intmR(r, mantaca, x_0, y_0, contactR)
sprintf("  expected num of points: E=%.2f", aaa)






# Plot average number of points against contact distribution as density func
contactRseq <- seq(from=1,to=7, length.out=50)
Fs <- array(dim=50)
Es <- array(dim=50)
numSim <- 40
k <- 1
for (contactR in contactRseq) {
  fallIn <- 0
  for (sim in seq(numSims)) {
    puntos <- rpoispp(mantaca, win = rectangulo)
    filled.contour(values, x = xtics, y = ytics, col=gray.colors(41),
                   plot.axes = { points(puntos); axis(1); axis(2) }
    )
    puntos_th <- rthin(X=puntos, P=rad_thin, r=r, nsim=1, drop=TRUE)
    filled.contour(values, x = xtics, y = ytics, col=gray.colors(41),
                   plot.axes = { points(puntos_th); axis(1); axis(2) }
    )
    # puntos_th <- puntos
    
    if (length(puntos_th$x) != 0) {
      # For the ocurrence of falling inside radius of contact
      for (i in seq(length(puntos_th$x))) {
        x <- puntos_th$x[i]
        y <- puntos_th$y[i]
        
        distance <- sqrt((x-x_0)^2 + (y-y_0)^2)
        if (distance <= contactR) {
          fallIn = fallIn + 1
          break
        }
      }
    }
  }
  
  Fs[k] <- fallIn / numSims
  Es[k] <- rad_th_intmR(r, mantaca, x_0, y_0, contactR)
  
  k <- k + 1
}

plot(contactRseq, Fs)
lines(contactRseq, Es)
