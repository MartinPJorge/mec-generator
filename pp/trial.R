library(spatstat)
library(graphics)

rm(list = ls())

gamma_lambda <- function(x, y, ... ){
  r <- sqrt(x^2 + y^2)
  return(dgamma(r, shape=5, scale=2))
}

# Define the rectangle region
xtics <- seq(-10, 10, length=10)
ytics <- seq(-10, 10, length=10)
values <- matrix(ncol = 10, nrow = 10)
i = 1;
for (x in xtics) {
  j = 1;
  for (y in ytics) {
    prob_ <- gamma_lambda(x, y)
    values[i, j] <- prob_
    j = j + 1
  }
  i = i + 1
}
rect_im <- im(values, xcol=xtics, yrow=ytics)

rectangulo <- owin(xrange=c(-10,10), yrange=c(-10,10))

puntos <- rpoispp(gamma_lambda, win = rectangulo)

# contour.im(rect_im, add=TRUE, axes=TRUE, clipwin = rectangulo, col=co)
#contour(rect_im, add=TRUE, axes=TRUE, clipwin = rectangulo)

filled.contour(values, x = xtics, y = ytics, col=gray.colors(40),
               plot.axes = { points(puntos); axis(1); axis(2) }
               )
#points(puntos$x, puntos$y)

# contour(rect_im, add=TRUE, axes=TRUE, clipwin = rectangulo)



