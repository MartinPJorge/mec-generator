library(ggmap)

madCenter <- c(left = -3.755549, bottom = 40.381710,
               right = -3.651810, top = 40.487886)
map <- get_map(madCenter, zoom = 13, maptype = "toner-lite")
ggmap(map)




# Approximation
e_ <- 8.1819190842622e-2 # first excentricity
a_ <- 6378137            # semi-major axis

madCenterR <- a_ * (1 - e_^2) / (1 - e_^2 * sin(40.44 *pi/180)^2)^(3/2)
madCenterW <- madCenterR * (-3.651810 - -3.755549)*pi/180*cos(40.381710 *pi/180)
