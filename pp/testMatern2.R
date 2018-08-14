library(spatstat)

rThinning <- 0.5

intensity <- function(x, y, m) {
  return(40 * y)
}


unitSquare <- owin(xrange = c(0, 10), yrange = c(0, 10))


points <- rMaternII(kappa = intensity, r = rThinning, win = unitSquare,
                    stationary=TRUE, nsim=1, drop=TRUE)
sprintf("> 5: %d", subset(points, y > 5)$n)
sprintf("< 5: %d", subset(points, y < 5)$n)
plot(x = points$x, y = points$y)
