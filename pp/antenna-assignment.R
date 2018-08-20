source("gen-utils-clean.R")

REGIONS <- "../data/regions.json"
REGIONS_FREQS <- "../data/frequencies.json"
REGION <- "Cobo-Calleja"


# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 1) {
    stop(cat("Arguments to receive are: ",
      "regionId\n"))
  } else {
    REGION <- args[1]
  }

# Get the grid JSONs
antennasGrids <- list.files(path = paste("../data/antennas/", REGION, sep = ""),
                           pattern = "gridded*")
peopleGrids <- list.files(path = paste("../data/people/", REGION, sep = ""),
                          pattern = "gridded*")
regionFreqs <- fromJSON(REGIONS_FREQS)

# Perform the assignment
for (a in 1:length(antennasGrids)) {
  aGrid <- fromJSON(antennasGrids[a])
  
  for (p in 1:length(peopleGrids)) {
    pGrid <- fromJSON(peopleGrids[p])
    if (pGrid$lonN == aGrid$lonN & pGrid$latN == aGrid$latN) {
      assignAntennas(antennasGrids[a], peopleGrids[p], regionFreqs)
    }
  }
}
