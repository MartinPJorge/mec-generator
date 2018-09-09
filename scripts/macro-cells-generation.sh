#!/bin/bash - 
#===============================================================================
#
#          FILE: macro-cells-generation.sh
# 
#         USAGE: ./macro-cells-generation.sh 
# 
#   DESCRIPTION: performs all the necessary steps to generate macro-cell
#                locations in an specific region.
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Jorge Martin Perez
#  ORGANIZATION: 
#       CREATED: 09/09/18 12:19
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

REGION_ID="Madrid-center"
LON_SAMPLES=100
LAT_SAMPLES=100
MACRO_CELL_FACTORS="../data/antennas/$REGION_ID/macro-cells-gen-factors.csv"
PEOPLE_INT_MAT="../data/people/$REGION_ID/people-intensity-matrix.csv"
PEOPLE_INT_LONS="../data/people/$REGION_ID/people-intensity-longitudes"
PEOPLE_INT_LATS="../data/people/$REGION_ID/people-intensity-latitudes"
MACRO_CELLS_DIR="../data/antennas/$REGION_ID/macro-cells-`date +%s`"

# Input script parameters
GENERATIONS=100 # number of macro-cells generations
AAUs=12 # number of AAUs per km^2



# Get inside the pp/ directory
cd ../pp

# Generate the people intensity matrix if needed
if [ ! -f $PEOPLE_INT_MAT ]; then
    echo "Running the people intensity matrix generation for region=$REGION_ID"
    Rscript --vanilla gen-people-lambda.R $REGION_ID $LON_SAMPLES $LAT_SAMPLES\
        $PEOPLE_INT_MAT $PEOPLE_INT_LONS $PEOPLE_INT_LATS
fi


# Generate the macro-cell factors if are not present
if [ ! -f $MACRO_CELL_FACTORS ]; then
    echo "Running the macro-cells factors generation for region=$REGION_ID"
    peopleMat longitudeSamp latitudeSamp AAUs outSquareFactors
    Rscript --vanilla macro-cells-factors.R $PEOPLE_INT_MAT $PEOPLE_INT_LONS\
       $PEOPLE_INT_LATS $AUUs $MACRO_CELL_FACTORS
fi


echo "Start the macro-cells generation"
mkdir -p $MACRO_CELLS_DIR
for i in `seq 1 $GENERATIONS`; do
    OUT_MACRO_CELLS="$MACRO_CELLS_DIR/macro-cells-$i.csv"
    echo "  macro-cell generation number $i"
    Rscript --vanilla macro-cells-generator.R $PEOPLE_INT_MAT $PEOPLE_INT_LONS\
       $PEOPLE_INT_LATS $MACRO_CELL_FACTORS $OUT_MACRO_CELLS
    echo "   output: $OUT_MACRO_CELLS"
done

