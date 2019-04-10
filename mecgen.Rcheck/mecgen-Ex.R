pkgname <- "mecgen"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "mecgen-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('mecgen')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cobo")
### * cobo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cobo
### Title: Generated antennas data for Cobo Calleja
### Aliases: cobo
### Keywords: datasets

### ** Examples

data(cobo)
assocs <- build5GScenario(lats = cobo$lat, lons = cobo$lon)

m1Assoc <- assocs[[1]]
m1Coords <- assocs[[2]]   
m1AccAssocs <- assocs[[3]]   
accCentCoords <- assocs[[4]]   
m2Assocs <- assocs[[5]]   
m2Switches <- assocs[[6]]   
m2AggAssocs <- assocs[[7]]   
aggCentCoords <- assocs[[8]]   
m3Assocs <- assocs[[9]]   
m3Switches <- assocs[[10]]  


write5GtoGraph(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                           m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                           m3Assocs, m3Switches, format = "gml",
               file = "/tmp/5g-mec.gml")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cobo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("regions")
### * regions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: regions
### Title: Regions data for the 5G infrastructure generation
### Aliases: regions
### Keywords: datasets

### ** Examples

data(regions)
regions$regions[[1]]$id
# [1] "Madrid-center"
regions$regions[[1]]$populations[[1]]$id
# [1] "Centro"
regions$regions[[1]]$populations[[1]]$population
# [1] 132284
regions$regions[[1]]$populations[[1]]$center
# [1] "40.416825,-3.706003" 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("regions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
