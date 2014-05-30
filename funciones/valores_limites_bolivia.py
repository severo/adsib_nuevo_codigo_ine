import os
from osgeo import ogr
import osr
from copy import deepcopy
import math
import pylab
import numpy as np

DATA_PATH = os.path.dirname("../datos/")
LIMITES_BOLIVIA_PATH = os.path.join(DATA_PATH, 'limite_nacional.gml')

# Abrir el archivo
ds = ogr.Open(LIMITES_BOLIVIA_PATH)
print ds.GetName()
print ds.GetLayerCount()

# Recuperar la primera y unica capa
layer = ds.GetLayer(0)
print layer.GetName()
print layer.GetFeatureCount()
origSpatialRef = layer.GetSpatialRef()

# Recuperar el unico "feature"
feature = layer.GetFeature(1)

# Definir la proyeccion conica conforme de Lambert de Bolivia
lambert = osr.SpatialReference()
lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Reproyectar en Lambert
SpatialRefTransform = osr.CoordinateTransformation(origSpatialRef,lambert)
geometry = feature.GetGeometryRef()
#print geometry.Length()
geometry.Transform(SpatialRefTransform)
#print geometry.Length()

# Calcular r y tita para cada punto
pts = geometry.GetGeometryRef(0)
points = []
x = []
y = []
r = []
tita = []
for p in xrange(pts.GetPointCount()):
	points.append((pts.GetX(p), pts.GetY(p)))
	x.append(pts.GetX(p))
	y.append(pts.GetY(p))
	points.append((pts.GetX(p), pts.GetY(p)))
	r.append(math.sqrt(pts.GetX(p)**2+ pts.GetY(p)**2))
	tita.append(np.arctan(pts.GetY(p)/pts.GetX(p)))

# Min y Max

print min(r)
print max(r)
print min(r)
print max(r)
print min(tita)
print max(tita)

print math.floor(min(r) * math.exp(min(tita)))
print math.floor(max(r) * math.exp(max(tita)))
