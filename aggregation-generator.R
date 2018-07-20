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
filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
               plot.axes = { points(puntos); axis(1); axis(2) }
)
filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
               plot.axes = { points(puntos_th); axis(1); axis(2) }
)
# points(antenas, pch=17)
# points(puntos, pch=19)
rectangulo_red <- owin(xrange=c(-10+r,10-r), yrange=c(-10+r,10-r))
int_th <- rad_th_intm(x=1, y=0, r=r, intensity=mantaca, obs_win=rectangulo)$Q
int_orig <- integral2(mantaca,
                 xmin = rectangulo$xrange[1], xmax = rectangulo$xrange[2],
                 ymin = rectangulo$yrange[1], ymax = rectangulo$yrange[2])$Q


cat("mu(R) = ", int_orig, "\n")
cat("mu_th(R) = ", int_th, "\n")
cat("mu(R) - mu_th(R) = ", int_orig - int_th, "\n")
cat("Realization X in R has ", length(puntos$x), " points\n")
cat("Realization X_th in R has ", length(puntos_th$x), " points\n")



