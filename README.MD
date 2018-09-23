# PPP scripts instructions
All the scripts files related to Poisson Point Processes are located under the `/pp` directory.

Every script must be executed as follows:
```bash
cd pp/
Rscript --vanilla script.R arg1 arg2 ...
```

Below it is shown a table summarizing the scripts' usage:

 name | description | usage
----- | ----------- | -----
along-rural-road.R | It generates a CSV with antennas along the A3 rural area road. The CSV with the rural area antennas must be generated. | `along-rural-road.R path/to/out.csv`
madrid-center-femto-cells.R | It generates a CSV with the femto-cells to be located in the Madrid skycraper | `madrid-center-femto-cells.R path/to/femto-cells.csv`
gen-people-lambda.R | It generates the people intensity function and store it as a matrix and two vectors that are the latitude and longitude samples. | `gen-people-lambda.R regionId longitudeSample latitudeSample matrixOUT lonsOUT latsOUT`
macro-cells-factors.R | Taking as input the people intensity matrix (and the corresponding longitude and latitude axis), it determines the parameters needed to obtain on average a certain amount of AAUs per km^2. See [macro-cells-factors note](#macro-cells-factors-generation) for output details. | `macro-cells-factors.R peopleMat longitudeSamp latitudeSamp AAUs outSquareFactors`
macro-cells-generator.R | It generates a CSV file with macro-cell antennas locations based in the associated factors per each km^2. It is based on the people intensity matrix and the factors obtained with `macro-cells-factors.R` to obtain the desired average AAUs with repulsion. | `macro-cells-generator.R peopleMat longitudeSamp latitudeSamp squareFactors macroCellsCSV`


## Script notes
### Macro-cells factors generation
* `factor / avgAAUs`: if people intensity function is multiplied by this at each km^2 square of the region, we get on average the AAUs passed by argument
* `sideSquares^2`: the amount of cells in which the km^2 is divided. The repulsion radius of generated macro-cells will be `1 / (2*sideSquares)`
* `smaller`: determines if the square area is 1km^2 or less. If set to TRUE, then it is a limiting rectangle of the region


# Automatic scripts
Under the `/script` directory there are some scripts that perform automatically specific tasks. Be sure you are inside the `/script` directory when you execute them.

script | description | usage
 ----- | ----------- | ----
`macro-cells-generation.sh` | It performs the necessary steps to generate macro-cells in the region contained in the `REGION_ID` variable inside the script. The generated macro-cells are located in a CSV under the `/data/antennas/REGION_ID/macro-cells-TIMESTAMP/` | `./macro-cells-generation.sh`
