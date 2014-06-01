import os
from osgeo import ogr
import osr
from copy import deepcopy
import math
import pylab
import numpy as np

# puntos de prueba - 1ra linea de la hoja de calculo "excel conversiones_codigo_unico.xlsx"
lat=-20.4491524218
lon=-64.3581217373

# sistema WGS84
wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)

# proyeccion conica conforme de Lambert de Bolivia
lambert = osr.SpatialReference()
lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# proyeccion WGFS84 -> Lambert
SpatialRefTransform = osr.CoordinateTransformation(wgs84,lambert)

# Convertir
x,y,z = SpatialRefTransform.TransformPoint(lon, lat, 0)

# Calcular r y tita
r = math.sqrt(x**2+ y**2)
tita = np.arctan(y/x) + np.pi

# Calcular z
z = math.floor(r * math.exp(tita))

# Min y Max

print "lat = "
print lat
print "lon = "
print lon
print "x = "
print x
print "y = "
print y
print "r = "
print r
print "tita = "
print tita
print "z = "
print z
