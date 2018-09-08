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
madrid-center-femto-cells.R | It generates a CSV with the femto-cells to be located in the Madrid skycraperswith the already existing LTE cells of operator net=7 | `las-femtooo.csv path/to/antennas-and-femto.csv`
gen-people-lambda.R | It generates the people intensity function and store it as a matrix and two vectors that are the latitude and longitude samples. | `gen-people-lambda.R regionId longitudeSample latitudeSample matrixOUT lonsOUT latsOUT`

