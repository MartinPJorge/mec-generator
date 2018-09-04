import csv
import json
import os

def add_i(row, i):
    if row[0:3] == range(row[0], row[0] +3 * i, i):
        return True

def getCoords(coordStr):
    """Obtains the coordinate longitude and latitude.

    :coordStr: coordinate string in format "12.3233,23.44343"
    :returns: (longitude, latitude)

    """
    coordsSpl = coordStr.split(",")
    return (float(coordsSpl[0]), float(coordsSpl[1]))


def fallInside(coordinates, region):
    """Tells if some coordinates fall inside a region.

    :coordinates: (longitude, latitude) as float numbers
    :region: { "id", "bl", "br", "tl", "tr" } with string coordinates
    :returns: Boolean

    """
    lat = coordinates[0]
    longi = coordinates[1]

    bl = getCoords(region["bl"])
    tr = getCoords(region["tr"])

    latitudeB = bl[0]
    longitudeL = bl[1]
    latitudeT = tr[0]
    longitudeR = tr[1]
    
    return longitudeL <= longi <= longitudeR and latitudeB <= lat <= latitudeT


def extractRegion(cellsCSV, region, fields):
    """Filters the cells falling inside the region from an incoming CSV file.

    :cellsCSV: CSV path to file containing cells from OpenCellID database.
    :region: { "id", "bl", "br", "tl", "tr" } with string coordinates
    :fields: list with the CSV fields for each column
    :returns: Nothing

    """
    regionId = str(region["id"]) 
    if not os.path.exists(regionId):
        os.mkdir(regionId)
    with open(regionId + "/" + regionId + ".csv", "w") as out_file:
        writer = csv.DictWriter(out_file, fields)
        writer.writeheader()
        with open(cellsCSV) as csv_file:
            for row in csv.DictReader(csv_file):
                rowCoords = [float(row["lat"]), float(row["lon"])]
                if fallInside(rowCoords, region):
                    writer.writerow(row)



CELLS_CSV = "spain-cells.csv"
REGIONS_FILE = "../regions.json"

if __name__ == "__main__":
    fields = []
    with open(CELLS_CSV) as f:
        fields = f.readline().split(",")
        fields[-1] = fields[-1][:-1] # remove \n

    with open(REGIONS_FILE) as f:
        regions = json.load(f)
        for region in regions["regions"]:
            extractRegion(CELLS_CSV, region, fields)

