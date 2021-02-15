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
mec-locations.R | It generates the locations of MEC PoPs based on the antennas coordinates and the radio technology used. Several methods for the selection of such locations are available. See [mec-locations note](#mec-locations) for more details. | `mec-locations.R regionID numLatSamples numLonSamples mecParamsJSON numMECs\<0 radioTech(FDD30kHz2sSPS, FDD120kHz7sSPS, TDD120kHz7sSPS) method(basicm1, basicm2, dig-outm1, dig-outm2, basicm1m2, dig-outm1m2) macro-cells.csv femto-cells.csv\NULL road-cells.csv\NULL MEClocations.csv intMatM1.csv\NULL intMatM2.csv\NULL`
plot-mec-locs.R | It plots the locations of the MEC PoPs with the M1 and M2 based intensity functions as heat-map along the locations. | `plot-mec-locs.R regionID latSamples lonSamples MEC_LOCATIONS_CSV MEC_M1_DIR\NULL MEC_M2_DIR\NULL OUT_PLOT_M1 OUT_PLOT_M2`

## Dependencies

The PPP scripts are written in [R](https://www.r-project.org/) and they use the following packages:

```
SDMTools
cluster
ggmap
graphics
latex2exp
metR
pracma
rjson
spatstat
stringi
```

## Script notes
### Macro-cells factors generation
* `factor / avgAAUs`: if people intensity function is multiplied by this at each km^2 square of the region, we get on average the AAUs passed by argument
* `sideSquares^2`: the amount of cells in which the km^2 is divided. The repulsion radius of generated macro-cells will be `1 / (2*sideSquares)`
* `smaller`: determines if the square area is 1km^2 or less. If set to TRUE, then it is a limiting rectangle of the region

### MEC locations
For the location of the MEC PoPs one can specify the NR radio technology used by he antennas:

carrier configuration | frequency | slot length | uplink scheduling | RTT
 -------------------- | --------- | ----------- | ----------------- | ---
 FDD | 30 kHz | 2 symbols | SPS | 0.78 ms
 FDD | 120 kHz | 7 symbols | SPS | 0.66 ms
 TDD | 120 kHz | 7 symbols | SPS | 0.78 ms
 
 **Note**: both first and last configurations give same results.
 
 As well several methods for the election of the MEC PoPs are available. The basic method is based on the reduction of the intensity function around the antennas once they are covered, whilst the dig-out is based on setting to zero the surroundings of the selected MEC PoP as well as the reduction around the covered antennas. The methods offer the possibility of selecting only locations at the M1, M2 ring or both. If both are selected, then it tries to minimize the number of MEC PoPs:

 *basicm1, basicm2, dig-outm1, dig-outm2, basicm1m2, dig-outm1m2*

# Automatic scripts
Under the `/script` directory there are some scripts that perform automatically specific tasks. Be sure you are inside the `/script` directory when you execute them.

script | description | usage
 ----- | ----------- | ----
`macro-cells-generation.sh` | It performs the necessary steps to generate macro-cells in the region contained in the `REGION_ID` variable inside the script. The generated macro-cells are located in a CSV under the `/data/antennas/REGION_ID/macro-cells-TIMESTAMP/` | `./macro-cells-generation.sh`



## Cite our work
If you use repo in your work, please consider citing our article:
```BibTeX
@ARTICLE{8667076,
  author={J. {Martín-Pérez} and L. {Cominardi} and C. J. {Bernardos} and A. {de la Oliva} and A. {Azcorra}},
  journal={IEEE Transactions on Broadcasting},
  title={Modeling Mobile Edge Computing Deployments for Low Latency Multimedia Services},
  year={2019},
  volume={65},
  number={2},
  pages={464-474},
  doi={10.1109/TBC.2019.2901406}
}
```

```txt
J. Martín-Pérez, L. Cominardi, C. J. Bernardos, A. de la Oliva and A. Azcorra, "Modeling Mobile Edge Computing Deployments for Low Latency Multimedia Services," in IEEE Transactions on Broadcasting, vol. 65, no. 2, pp. 464-474, June 2019, doi: 10.1109/TBC.2019.2901406.
```
