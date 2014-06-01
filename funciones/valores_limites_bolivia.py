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
print "Archivo:"
print ds.GetName()

# Recuperar la primera y unica capa
layer = ds.GetLayer(0)
print "Capa: "
print layer.GetName()
origSpatialRef = layer.GetSpatialRef()

# Recuperar la geometria del unico "feature"
feature = layer.GetFeature(1)
geometry = feature.GetGeometryRef()

# Definir la proyeccion conica conforme de Lambert de Bolivia
lambert = osr.SpatialReference()
lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Reproyectar en Lambert
SpatialRefTransform = osr.CoordinateTransformation(origSpatialRef,lambert)
geometry.Transform(SpatialRefTransform)

# Calcular r y tita para cada punto
pts = geometry.GetGeometryRef(0)
vx = []
vy = []
vr = []
vtita = []
vz = []
for p in xrange(pts.GetPointCount()):
	x = pts.GetX(p)
	y = pts.GetY(p)
	r = math.sqrt(x**2+ y**2)
	tita = np.arctan(y/x)+np.pi
	z = math.floor(r * math.exp(tita))
	
	vx.append(x)
	vy.append(y)
	vr.append(r)
	vtita.append(tita)
	vz.append(z)

# Min y Max

print "Valores de x"
print min(vx)
print max(vx)
print "Valores de y"
print min(vy)
print max(vy)
print "Valores de r"
print min(vr)
print max(vr)
print "Valores de tita"
print min(vtita)
print max(vtita)
print "Valores de z"
print min(vz)
print max(vz)
#print math.floor(min(vr) * math.exp(min(vtita)))
#print math.floor(max(vr) * math.exp(max(vtita)))

# Numero de codigos
print "Numero de codigos"
print max(vz) - min(vz)

# Superficie de Bolivia
area_bolivia = geometry.ConvexHull().GetArea()
print "Superficie Bolivia (m2) = %d" % area_bolivia
