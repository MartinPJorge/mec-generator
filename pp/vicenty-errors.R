library('SDMTools')

madridBl = c(40.487886, -3.755549)
madridTr = c(40.487886, -3.651810)
dirAngle = 90

# Point in the middle of the region
madridCent = c(madridTr[1] - madridBl[1], madridTr[2] - madridBl[2])
madridCent = c(madridBl[1] + madridCent[1] / 2, madridBl[2] + madridCent[2] / 2)

# Point to the right in 500 meters
endP = destination(lat=madridCent[1], lon=madridCent[2], bearing=90, distance = 500)

# 500m in longitude is v degree
longiDeg = endP$lon2 - madridCent[2]
distance(lat1=madridCent[1], lon1=madridCent[2],
         lat2=madridCent[1]+longiDeg, lon2=madridCent[2])

