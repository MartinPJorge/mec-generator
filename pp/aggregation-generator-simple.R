rm(list = ls())
source("gen-utils-simple.R")
library(spatstat)

library(graphics) 
  
graphics.off()


intensity_f <- function(x, y) {
  return(y)
}


xtics <- seq(0, 1, length=100)
ytics <- seq(0, 1, length=100)
values <- matrix(ncol=100, nrow=100)
i = 1;
for (x in xtics) {
  j = 1;
  for (y in ytics) {
    values[i, j] <- intensity_f(x, y)
    j = j + 1
  }
  i = i + 1
}


r = 0.1 # thinning radius
rectangulo <- owin(xrange=c(0,1), yrange=c(0,1))
puntos <- rpoispp(intensity_f, win = rectangulo)
#contour.im(rect_im, add=FALSE, axes=TRUE, clipwin = rectangulo)
#contour(rect_im, add=TRUE, axes=TRUE, clipwin = rectangulo)
puntos_th_own <- rthin(X=puntos, P=rad_thin, r=r, nsim=1, drop=TRUE)
puntos_th <- rMaternI(kappa = intensity_f, r=r, win = rectangulo)

abs(length(puntos_th_own$x) - length(puntos_th$x))

filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
               plot.axes = { points(puntos); axis(1); axis(2) }
)
filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
              plot.axes = { points(puntos_th); axis(1); axis(2) }
)



# points(antenas, pch=17)
# points(puntos, pch=19)
rectangulo <- owin(xrange=c(0.4, 0.6), yrange=c(0.2,0.4))
int_th <- rad_th_intm(x=1, y=0, r=r, intensity=intensity_f, obs_win=rectangulo)$Q
int_orig <- integral2(intensity_f,
                 xmin = rectangulo$xrange[1], xmax = rectangulo$xrange[2],
                 ymin = rectangulo$yrange[1], ymax = rectangulo$yrange[2])$Q


cat("mu(R) = ", int_orig, "\n")
cat("mu_th(R) = ", int_th, "\n")
cat("mu(R) - mu_th(R) = ", int_orig - int_th, "\n")
cat("Realization X in R has ", length(puntos$x), " points\n")
cat("Realization X_th in R has ", length(puntos_th$x), " points\n")








########## CHECK
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

ball_intm(0,1, 1, intensity_f)