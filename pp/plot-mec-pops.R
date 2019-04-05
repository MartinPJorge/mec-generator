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
library(stringi)

########### MADRID PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Madrid-center"
MEC_LOCATIONS_CSV <- paste("../data/mec-pops/Madrid-center/mec-locations-i-iii",
                           "/Madrid-Center-12AAUs-factor16-1.csv", sep = "")
MEC_M1_DIR <- "../data/mec-pops/Madrid-center/mec-m1-mats-i-iii"
MEC_M2_DIR <- "../data/mec-pops/Madrid-center/mec-m2-mats-i-iii"
PEOPLE_LONS <- "../data/people/Madrid-centro/people-intensity-longitudes"
PEOPLE_LATS <- "../data/people/Madrid-centro/people-intensity-latitudes"
CLI <- FALSE # flag to tell if file is executed from CLI
OUT_PLOT_M1 <- NULL 
OUT_PLOT_M2 <- NULL
#####################################


########### COBO CALLEJA PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Cobo-Calleja"
MEC_LOCATIONS_CSV <- "../data/mec-pops/Cobo-Calleja/mec-locations-i-iii/Cobo-Calleja-12AAUs-factor1.csv"
MEC_M1_DIR <- "../data/mec-pops/Cobo-Calleja/mec-m1-mats-i-iii"
MEC_M2_DIR <- "../data/mec-pops/Cobo-Calleja/mec-m2-mats-i-iii"
PEOPLE_LONS <- "../data/people/Cobo-Calleja/people-intensity-lons"
PEOPLE_LATS <- "../data/people/Cobo-Calleja/people-intensity-lats"
CLI <- FALSE # flag to tell if file is executed from CLI
OUT_PLOT_M1 <- NULL 
OUT_PLOT_M2 <- NULL
###########################################



########### HOCES DEL CABRIEL PARAMS ###########
REGIONS <- "../data/regions.json"
REGION_NAME <- "Hoces-del-Cabriel"
MEC_LOCATIONS_CSV <- paste("../data/mec-pops/Hoces-del-Cabriel/",
                           "mec-locations-i-iii-road/",
                           "Hoces-del-Cabriel-road-1.csv", sep = "")
MEC_M1_DIR <- "../data/mec-pops/Hoces-del-Cabriel/mec-m1-mats-i-iii-road"
MEC_M2_DIR <- "../data/mec-pops/Hoces-del-Cabriel/mec-m2-mats-i-iii-road"
PEOPLE_LONS <- "../data/people/Hoces-del-Cabriel/people-intensity-lons"
PEOPLE_LATS <- "../data/people/Hoces-del-Cabriel/people-intensity-lats"
CLI <- FALSE # flag to tell if file is executed from CLI
OUT_PLOT_M1 <- NULL 
OUT_PLOT_M2 <- NULL
################################################



# Get the radio technology
radioTech = ifelse(stri_detect(MEC_M1_DIR, regex="iii"),
                   yes = "FDD 30kHz 2 symbols SPS, TDD 120kHz 7 symbols SPS",
                   no = "FDD 120kHz 7 symbols SPS")


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
mecLocs$circleSize <- 1 + ceil(10 * mecLocs$coveredAs / max(mecLocs$coveredAs))
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



# Get the subregion to be plotted
if (!is.null(region$plotDetails$tr)) {
  region$bl <- getCoords(region$plotDetails$bl)
  region$br <- getCoords(region$plotDetails$br)
  region$tl <- getCoords(region$plotDetails$tl)
  region$tr <- getCoords(region$plotDetails$tr)
  mecLocs <- subset(mecLocs, mecLocs$lon > region$bl$lon &
                      mecLocs$lon < region$br$lon &
                      mecLocs$lat > region$br$lat & mecLocs$lat < region$tr$lat)
}


# Obtain the average of all the M1 matrices
avgM1 <- NULL
m1MatLons <- lonAxis
m1MatLats <- latAxis
if (!is.null(MEC_M1_DIR)) {
  m1s <- 0
  
  for (M1_CSV in list.files(path = MEC_M1_DIR, full.names = TRUE)) {
    chop <- NULL
    readM1 <- as.matrix(read.csv(file = M1_CSV)[,-1])
      
    if (!is.null(region$plotDetails$tr)) { # chopped region
      chop <- chopIntMat(intMat = readM1, lonAxis = lonAxis, latAxis = latAxis,
                         lonL = region$bl$lon, lonR = region$tr$lon,
                         latB = region$bl$lat, latT = region$tr$lat)
      if (is.null(avgM1)) {
        avgM1 <- chop$mat
        m1MatLats <- chop$latAxis
        m1MatLons <- chop$lonAxis
      } else {
        avgM1 <- avgM1 + chop$mat
      }
    } else { # full region
      if (is.null(avgM1)) {
        avgM1 <- matrix(data = 0, nrow = length(lonAxis),
                        ncol = length(latAxis))
      }
      avgM1 <- avgM1 + readM1
    }
    
    m1s <- m1s + 1
  }
  
  avgM1 <- avgM1 / m1s
}


# Obtain the average of all the M2 matrices
avgM2 <- NULL
m2MatLons <- lonAxis
m2MatLats <- latAxis
if (!is.null(MEC_M2_DIR)) {
  m2s <- 0
  
  for (M2_CSV in list.files(path = MEC_M2_DIR, full.names = TRUE)) {
    chop <- NULL
    readM2 <- as.matrix(read.csv(file = M2_CSV)[,-1])
      
    if (!is.null(region$plotDetails$tr)) { # chopped region
      chop <- chopIntMat(intMat = readM2, lonAxis = lonAxis, latAxis = latAxis,
                         lonL = region$bl$lon, lonR = region$tr$lon,
                         latB = region$bl$lat, latT = region$tr$lat)
      if (is.null(avgM2)) {
        avgM2 <- chop$mat
        m2MatLats <- chop$latAxis
        m2MatLons <- chop$lonAxis
      } else {
        avgM2 <- avgM2 + chop$mat
      }
    } else { # full region
      if (is.null(avgM2)) {
        avgM2 <- matrix(data = 0, nrow = length(lonAxis),
                        ncol = length(latAxis))
      }
      avgM2 <- avgM2 + readM2
    }
    
    m2s <- m2s + 1
  }
  
  avgM2 <- avgM2 / m2s
}




# Get the map
mapRegion <- c(left = region$bl$lon, bottom = region$bl$lat,
               right = region$tr$lon, top = region$tr$lat)
map <- get_map(mapRegion, zoom = region$plotDetails$zoom, source = "stamen",
              maptype = region$plotDetails$mapType)

# Plot
if (!is.null(MEC_M1_DIR)) {
  # Get the groups correctly for the contour
  levels <- pretty(range(avgM1), n = 7)
  p <- contourPolys::fcontour(x = m1MatLons, y = m1MatLats,
                              z = avgM1, levels)
  m <- cbind(x = unlist(p[[1]]), 
             y = unlist(p[[2]]), 
             lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
             upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
             g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 
  
  gd <- as.data.frame(m)
  mynamestheme <- theme(plot.title = element_text(family = "Helvetica",
                                                  face = "bold", size = (28)),
                        axis.title = element_text(size = 18),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18),
                        legend.text = element_text(size = 20),
                        legend.title = element_text(size = 28)
                        )
  
  
  # Plot
  gg_m <- ggmap(map)
  gg_m <- gg_m + geom_polygon(data=gd, aes(x, y, group = g, fill = upper, alpha = upper)) +
    scale_fill_gradient(low = "gray", high = "black") +
    scale_alpha(range = c(0, 0.4)) +
    geom_point(data=mecLocs, aes(x=lon, y=lat, shape=shapes,
                                     size=circleSize)) +
    scale_shape_manual(name = "Network ring",
                       labels = c("M1", "M2"),
                       values = c(18, 15)) +
    labs(size = TeX("$\\frac{10·C_A(mp)}{max_{mp \\in PoPs} C_A(mp)}$")) +
    labs(alpha = TeX("C_{A,M1}"), x = "LONGITUDE", y = "LATITUDE") +
    guides(fill = FALSE, size = FALSE) +
    mynamestheme
    #ggtitle(label = "MEC PoPs", subtitle = sprintf("antennas' radio: %s", radioTech)) +
    gg_m
  
  if (CLI) {
    ggsave(filename = OUT_PLOT_M1, plot = gg_m)
  }
}
if (!is.null(MEC_M2_DIR)) {
  # Get the groups correctly for the contour
  levels <- pretty(range(avgM2), n = 7)
  p <- contourPolys::fcontour(x = m2MatLons, y = m2MatLats,
                              z = avgM2, levels)
  m <- cbind(x = unlist(p[[1]]), 
             y = unlist(p[[2]]), 
             lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
             upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
             g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 
  
  gd <- as.data.frame(m)
  mynamestheme <- theme(plot.title = element_text(family = "Helvetica",
                                                  face = "bold", size = (28)),
                        axis.title = element_text(size = 18),
                        axis.text.x = element_text(size = 18),
                        axis.text.y = element_text(size = 18),
                        legend.text = element_text(size = 20),
                        legend.title = element_text(size = 28)
                        )
  
  
  # Plot
  gg_m <- ggmap(map)
  gg_m <- gg_m + geom_polygon(data=gd, aes(x, y, group = g, fill = upper, alpha = upper)) +
    scale_fill_gradient(low = "gray", high = "black") +
    scale_alpha(range = c(0.3, 0.6)) +
    geom_point(data=mecLocs, aes(x=lon, y=lat, shape=shapes,
                                     size=circleSize)) +
    labs(size = TeX("$\\frac{10·C_A(mp)}{max_{mp \\in PoPs} C_A(mp)}$")) +
    labs(alpha = TeX("C_{A,M2}"), x = "LONGITUDE", y = "LATITUDE") +
    guides(fill = FALSE, size = FALSE) +
    scale_shape_manual(name = "Network ring",
                       labels = c("M2", "M1"),
                       values = c(15, 18)) +
    mynamestheme
    #ggtitle(label = "MEC PoPs", subtitle = sprintf("antennas' radio: %s", radioTech)) +
    gg_m
  
    #ggtitle("MEC PoPs", subtitle = sprintf("antennas' radio: %s", radioTech)) +
    #ggtitle(label = "MEC PoPs", subtitle = sprintf("antennas' radio: %s", radioTech)) +
  
  if (CLI) {
    ggsave(filename = OUT_PLOT_M2, plot = gg_m)
  }
}

