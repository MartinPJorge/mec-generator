# mec-generator
This branch is a R package `mecgenerator`.
Its purpose is the automatic generation of 5G MEC infrastructure graphs
as the ones depicted in figure below [1], given a map region:
![Image of the depicted 5G MEC infrastructure with antennas, and up to the core layer](https://github.com/MartinPJorge/mec-generator/blob/5g-infra-gen/img/infra.gif)

To install it just clone the repo under the desired location and run inside R:
```R
setwd("/path/to/mec-generator")
devtools::install("mecgen")
```
# Usage
## Cells generation
TODO - add reference to article methodology

## Rings generation
The package comes with a `data.frame` with antennas' locations inside the area of Cobo Calleja, an industrial area in the
south of Madrid.
Given the set of antennas, it generates the whole switches of the shown figure, along with the switches and interconnections among them.
The following code:
```R
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
```

generates graphs as the shown below, which conveys the architecture of first picture:
![Image with the graph representation of the 5G MEC infrastructure](https://github.com/MartinPJorge/mec-generator/blob/5g-infra-gen/img/infra-graphs.png)


### References
[1] L. Cominardi, L. M. Contreras, C. J. Bcrnardos and I. Berberana, **"Understanding QoS Applicability in 5G Transport Networks,"** 2018 IEEE International Symposium on Broadband Multimedia Systems and Broadcasting (BMSB), Valencia, 2018, pp. 1-5.
doi: 10.1109/BMSB.2018.8436847
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8436847&isnumber=8436592

[2] J. Martín-Pérez, L. Cominardi, C. J. Bernardos, A. de la Oliva and A. Azcorra, **"Modeling Mobile Edge Computing Deployments for Low Latency Multimedia Services,"** in IEEE Transactions on Broadcasting.
doi: 10.1109/TBC.2019.2901406
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8667076&isnumber=4359969
