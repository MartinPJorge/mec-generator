source("gen-utils-clean.R")
library(pracma)
library(SDMTools)
library(ggmap)
library(spatstat)
library(graphics) 
library(rjson)
library(cluster)
library(latex2exp)
library(metR)

########### MADRID PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
MEC_LOCATIONS_CSV <- paste("../data/mec-pops/Madrid-center/mec-locations-ii",
                           "/Madrid-center-12AAUs-factor16-1.csv", sep = "")
MEC_M1_DIR <- "../data/mec-pops/Madrid-center/mec-m1-mats-ii"
MEC_M2_DIR <- "../data/mec-pops/Madrid-center/mec-m2-mats-ii"
PEOPLE_LONS <- "../data/people/Madrid-centro/people-intensity-longitudes"
PEOPLE_LATS <- "../data/people/Madrid-centro/people-intensity-latitudes"
CLI <- FALSE # flag to tell if file is executed from CLI
OUT_PLOT_M1 <- NULL 
OUT_PLOT_M2 <- NULL
#####################################


########### COBO CALLEJA PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Cobo-Calleja"
MEC_LOCATIONS_CSV <- "../data/mec-pops/Cobo-Calleja/mec-locations-i-iii/Cobo-Calleja-12AAUs-factor13-pops-1.csv"
MEC_M1_DIR <- "../data/mec-pops/Cobo-Calleja/mec-m1-mats-i-iii"
MEC_M2_DIR <- "../data/mec-pops/Cobo-Calleja/mec-m2-mats-i-iii"
PEOPLE_LONS <- "../data/people/Cobo-Calleja/people-intensity-lons"
PEOPLE_LATS <- "../data/people/Cobo-Calleja/people-intensity-lats"
CLI <- FALSE # flag to tell if file is executed from CLI
OUT_PLOT_M1 <- NULL 
OUT_PLOT_M2 <- NULL
###########################################




# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 8) {
    stop(paste("Arguments to receive are: ",
      "regionID latSamples lonSamples MEC_LOCATIONS_CSV MEC_M1_DIR|NULL",
      "MEC_M2_DIR|NULL OUT_PLOT_M1 OUT_PLOT_M2"
    ))
  } else {
    CLI <- TRUE
    REGION_NAME <- args[1]
    PEOPLE_LATS <- args[2]
    PEOPLE_LONS <- args[3]
    MEC_LOCATIONS_CSV <- args[4]
    MEC_M1_DIR <- args[5]
    if (MEC_M1_DIR == "NULL") {
      MEC_M1_DIR = NULL
    }
    MEC_M2_DIR <- args[6]
    if (MEC_M2_DIR == "NULL") {
      MEC_M2_DIR = NULL
    }
    OUT_PLOT_M1 <- args[7]
    OUT_PLOT_M2 <- args[8]
    
    if (is.null(MEC_M1_DIR) & is.null(MEC_M2_DIR)) {
      stop("M1 or M2 matrices must be provided, but none is given")
    }
  }


# Load files
regions <- fromJSON(file = REGIONS)
lonAxis <- scan(file = PEOPLE_LONS)
latAxis <- scan(file = PEOPLE_LATS)
mecLocs <- read.csv(file = MEC_LOCATIONS_CSV)
mecLocs$circleSize <- ceil(10 * mecLocs$coveredAs / max(mecLocs$coveredAs))
shapes <- c()
for (row in 1:nrow(mecLocs)) {
  shapes <- c(shapes, ifelse(mecLocs[row,]$ring == "M2", 15, 18))
}
mecLocs$shapes <- as.factor(shapes)

# Find the selected region
region <- NULL
for (region_ in regions$regions) {
  if (region_$id == REGION_NAME) {
    region <- list(id = region_$id,
                bl = getCoords(region_$bl),
                br = getCoords(region_$br),
                tl = getCoords(region_$tl),
                tr = getCoords(region_$tr),
                repulsionRadius = region_$repulsionRadius,
                plotDetails = region_$plotDetails,
                populations = region_$populations)
    break
  }
}


# Obtain the average of all the M1 matrices
avgM1 <- NULL
avgM1Fr <- NULL
if (!is.null(MEC_M1_DIR)) {
  avgM1 <- matrix(data = 0, nrow = length(lonAxis), ncol = length(latAxis))
  
  m1s <- 0
  for (M1_CSV in list.files(path = MEC_M1_DIR, full.names = TRUE)) {
    avgM1 <- avgM1 + as.matrix(read.csv(file = M1_CSV)[,-1])
    m1s <- m1s + 1
  }
  
  avgM1 <- avgM1 / m1s
  avgM1Fr <- intMat2Frame(intMat = avgM1,
                           latAxis = latAxis,
                           lonAxis = lonAxis)
}


# Obtain the average of all the M1 matrices
avgM2 <- NULL
avgM2Fr <- NULL
if (!is.null(MEC_M2_DIR)) {
  avgM2 <- matrix(data = 0, nrow = length(lonAxis), ncol = length(latAxis))
  
  m2s <- 0
  for (M2_CSV in list.files(path = MEC_M2_DIR, full.names = TRUE)) {
    avgM2 <- avgM2 + as.matrix(read.csv(file = M2_CSV)[,-1])
    m2s <- m2s + 2
  }
  
  avgM2 <- avgM2 / m2s
  avgM2Fr <- intMat2Frame(intMat = avgM2,
                           latAxis = latAxis,
                           lonAxis = lonAxis)
}



# Get the map
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
              maptype = region$plotDetails$mapType)

# Plot
if (!is.null(MEC_M1_DIR)) {
  # Get the groups correctly for the contour
  levels <- pretty(range(avgM1Fr$intensity), n = 7)
  p <- contourPolys::fcontour(x = lonAxis, y = latAxis,
                              z = avgM1, levels)
  m <- cbind(x = unlist(p[[1]]), 
             y = unlist(p[[2]]), 
             lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
             upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
             g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 
  
  gd <- as.data.frame(m)
  
  # Plot
  gg_m <- ggmap(map)
  gg_m <- gg_m + geom_polygon(data=gd, aes(x, y, group = g, fill = upper, alpha = upper)) +
    scale_fill_gradient(low = "gray", high = "black") +
    scale_alpha(range = c(0.3, 0.6)) +
    geom_point(data=mecLocs, aes(x=lon, y=lat, shape=shapes,
                                     size=circleSize)) +
    labs(size = TeX("$\\frac{10·C_A(mp)}{max_{mp \\in PoPs} C_A(mp)}$")) +
    labs(alpha = TeX("C_A(mp)"), x = "longitude", y = "latitude") +
    guides(fill = FALSE) +
    scale_shape_manual(name = "Network ring",
                       labels = c("M2", "M1"),
                       values = c(15, 18))
  gg_m
  
  if (CLI) {
    ggsave(filename = OUT_PLOT_M1, plot = gg_m)
  }
}
if (!is.null(MEC_M2_DIR)) {
  # Get the groups correctly for the contour
  levels <- pretty(range(avgM2Fr$intensity), n = 7)
  p <- contourPolys::fcontour(x = lonAxis, y = latAxis,
                              z = avgM2, levels)
  m <- cbind(x = unlist(p[[1]]), 
             y = unlist(p[[2]]), 
             lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
             upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
             g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 
  
  gd <- as.data.frame(m)
  
  # Plot
  gg_m <- ggmap(map)
  gg_m <- gg_m + geom_polygon(data=gd, aes(x, y, group = g, fill = upper, alpha = upper)) +
    scale_fill_gradient(low = "grey45", high = "black") +
    scale_alpha(range = c(0.3, 0.9)) +
    geom_point(data=mecLocs, aes(x=lon, y=lat, shape=shapes,
                                     size=circleSize)) +
    labs(size = TeX("$\\frac{10·C_A(mp)}{max_{mp \\in PoPs} C_A(mp)}$")) +
    labs(alpha = TeX("C_A(mp)"), x = "longitude", y = "latitude") +
    guides(fill = FALSE) +
    scale_shape_manual(name = "Network ring",
                       labels = c("M2", "M1"),
                       values = c(15, 18))
  gg_m
  
  if (CLI) {
    ggsave(filename = OUT_PLOT_M2, plot = gg_m)
  }
}

