 mec-generator
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
The package contains a `data.frame` with antennas' locations inside the area of Cobo Calleja, an industrial area in the
south of Madrid.
Given the set of antennas, it generates the whole switches of the shown figure, along with the switches and interconnections among them.
The following code:
```R
library(mecgen)

coboCells <- mecgen::cobo
regions <- mecgen::regions
assocs <- build5GScenario(lats = coboCells$lat, lons = coboCells$lon)
```
creates a list of the followinig `data.frames`:



|   list item  | variable      | meaning                                                                            |
|:------------:|---------------|------------------------------------------------------------------------------------|
| assocs[[1]]  | m1Assoc       | Cell antenna coordinates and M1 switch id.                                         |
| assocs[[2]]  | m1Coords      | M1 switches coordinates and their ids.                                             |
| assocs[[3]]  | m1AccAssocs   | M1 coordinates and their access ring ids.                                          |
| assocs[[4]]  | accCentCoords | Coordinates of the access ring centers.                                            |
| assocs[[5]]  | m2Assocs      | Associations of access ring centers to the M2 switches' ids.                       |
| assocs[[6]]  | m2Switches    | Coordinates of M2 switches and their ids.                                          |
| assocs[[7]]  | m2AggAssocs   | M2 switches coordinates and ids of aggregation ring to which they are associated.  |
| assocs[[8]]  | aggCentCoords | Coordinates of aggregation rings with their ids.                                   |
| assocs[[9]]  | m3Assocs      | Coordinates of aggregation ring centers, and the M3 switch they are associated to. |
| assocs[[10]] | m3Switches    | Coordinates of M3 switches with their ids.                                         |


These `data.frame`s are used to create another two `data.frame`s containing the nodes and edges of the infrastructure graph :


```R
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


frames <- graphFrames(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                           m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                           m3Assocs, m3Switches)
```
the result contains `frames$nodes` and `frames$links`. Both of them are `data.frame`s with the infrastructure nodes and interconnections. Both of them can be latter on modified to include servers.

## Generate compute nodes and endpoints

### Server generation
To generate servers attached to M1, M2 or M3 switches, just execute `attachServers()` function with the retrieved `data.frame`s of the infrastructure. The following chunk of code generates 3 servers with 12 Mbps connections to M2 switches, and locates them 0 meters away:

```R
attachFrames <- attachServers(nodes = frames$nodes, links = frames$links,
                              numServers = 3,
                              bandwidth = 12,
                              bandwidthUnits = "Mbps",
                              distance = 0,
                              distanceUnits = "meter",
                              switchType = "m2",
                              properties = list(cpu=2, mem=20, disk=100), idPrefix = "dell")
```
as well it provides a list of properties that each server has. The `properties` parameter can be any list of characteristics that generated servers will have.

**Note**: the `idPrefix` parameter must be unique for every execution of `attachServer()`, otherwise the library can't distinguish new servers from previous ones with same prefix.

### Fog nodes generation
For the generation of fog compute nodes, just execute  `attachFogNodes` as `attachServers()`, but providing the limiting coordinates for the coordinates of the fog nodes that are generated. The function will attach the fog nodes to the closes cell antenna inside the `nodes` `data.frame` parameter.

```R
# Retrieve the Cobo Calleja region limiting coordinates
coboRegion <- regions$regions[[2]]
coboLonL <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][2])
coboLonR <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][2])
coboLatB <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][1])
coboLatT <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][1])

attachFrames <- attachFogNodes(nodes = attachFrames$nodes,
                               links = attachFrames$links,
                               latB = coboLatB, latT = coboLatT,
                               lonL = coboLonL, lonR = coboLonR,
                               numNodes = 10,
                               properties = list(cpu = 2, mem = 20), bandwidth = 20,
                               bandwidthUnits = "Mpbs",
                               idPrefix = "test")
```

### Endpoints generation
Endpoints can represent users consuming infrastructure resources, and `attachEndpoints()` generates them and associates them to the nearest cell antennas.
```R
attachFrames <- attachEndpoints(nodes = attachFrames$nodes,
                                links = attachFrames$links, latB = coboLatB,
                                latT = coboLatT, lonL = coboLonL,
                                lonR = coboLonR, numEndpoints = 10,
                                bandwidth = 10,
                                bandwidthUnits = "Mbps",
                                idPrefix = "person")

```

## Graph creation
The nodes and links `data.frame`s present in the result of `graphFrames()`, `attachServers()`, `attachFogNodes()`, and `attachEndpoints()`; are used as input for the generation of undirected graphs using `igraph` package:
```R
links <- attachFrames$links
nodes <- attachFrames$nodes

g = igraph::graph_from_data_frame(links, vertices = nodes, directed = FALSE)
igraph::write_graph(graph = g, file = "/tmp/5g-mec.gml", format = "gml")
```


Picture below represents the graph representation of the 5G MEC infrastructure (without server, fog nodes, and enpoint attachments):
![Image with the graph representation of the 5G MEC infrastructure](https://github.com/MartinPJorge/mec-generator/blob/5g-infra-gen/img/infra-graphs.png)


### References
[1] L. Cominardi, L. M. Contreras, C. J. Bcrnardos and I. Berberana, **"Understanding QoS Applicability in 5G Transport Networks,"** 2018 IEEE International Symposium on Broadband Multimedia Systems and Broadcasting (BMSB), Valencia, 2018, pp. 1-5.
doi: 10.1109/BMSB.2018.8436847
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8436847&isnumber=8436592

[2] J. Martín-Pérez, L. Cominardi, C. J. Bernardos, A. de la Oliva and A. Azcorra, **"Modeling Mobile Edge Computing Deployments for Low Latency Multimedia Services,"** in IEEE Transactions on Broadcasting.
doi: 10.1109/TBC.2019.2901406
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8667076&isnumber=4359969


