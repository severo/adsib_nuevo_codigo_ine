# Para hacer funcionar con Debian
#  sudo aptitude install python-geohash python-gdal python-pip
# instalar la libreria python mgrs:
#  sudo pip install mgrs

import os
from osgeo import ogr
import osr
from copy import deepcopy
import geohash
import math
import pylab
import numpy as np

def calcular_codigo(lon,lat):
	# Definir la proyeccion conica conforme de Lambert de Bolivia
	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

	# Reproyectar en Lambert
	SpatialRefTransform = osr.CoordinateTransformation(origSpatialRef,lambert)
	x,y,z = SpatialRefTransform.TransformPoint(lon, lat, 0)

	# Calcular r y tita
	r = math.sqrt(x**2+ y**2)
	tita = np.arctan(y/x)+np.pi
	z = math.floor(r * math.exp(tita))

	return (r,tita,z)

# Calcula el codigo INE que pasa por un punto
def calcular_codigo_ine(lon,lat,nombre,poly_bolivia):
	(r,tita,z) = calcular_codigo(lon,lat)

	## Generar un SHP con esta franja
	drv = ogr.GetDriverByName("ESRI Shapefile")
	outputfile = "../resultados/INE_" + nombre + ".shp"
	if os.path.exists(outputfile):
		drv.DeleteDataSource(outputfile)
	output = drv.CreateDataSource(outputfile)
	franja = output.CreateLayer("Franja", geom_type=ogr.wkbPolygon)

	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
	wgs84 = osr.SpatialReference()
	wgs84.ImportFromEPSG(4326)
	transfWGS84 = osr.CoordinateTransformation(lambert,wgs84)

	## Calcular la franja de territorio con el mismo codigo que el punto
	mintita = min(vtita)
	maxtita = max(vtita)
	paso_tita = (maxtita-mintita)/1000
	vrmin = []
	vptmin = []
	vrmax = []
	vptmax = []
	vdeltar = []
	for i in range(0,999):
		tita = mintita + i*paso_tita
		rmin = z * math.exp(- tita)
		rmax = (z+1) * math.exp(- tita)
		vrmin.append(rmin)
		vptmin.append((rmin,tita))
		vrmax.append(rmax)
		vptmax.append((rmax,tita))
		vdeltar.append(rmax-rmin)

	vpt = vptmin
	vpt.extend(reversed(vptmax))

	ring = ogr.Geometry(ogr.wkbLinearRing)
	for i in vpt:
		x_franja = i[0]*math.cos(i[1]-math.pi)
		y_franja = i[0]*math.sin(i[1]-math.pi)
		ring.AddPoint(x_franja,y_franja)

	ring.CloseRings()

	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)

	poly_intersection = poly.Intersection(poly_bolivia)

	print "************************************************"
	print "Codigo INE z - " + nombre + ":"
	print " * valor: %d" % z
	print " * area equi-codigo:"
	print "   * superficie: %fm2 (equivalente de un cuadrado de %fm de lado)" % (poly_intersection.GetArea(),math.sqrt(poly_intersection.GetArea()))
	print "   * es un 'arco' de %fkm de largo" % (poly_intersection.Boundary().Length()/2)
	print "                   y %fm de ancho (promedio - min=%fm max=%fm)" % (float( sum(vdeltar) ) / len(vdeltar),min(vdeltar),max(vdeltar))

	poly_intersection.Transform(transfWGS84)
	f = ogr.Feature(feature_def=franja.GetLayerDefn())
	f.SetGeometryDirectly(poly_intersection)
	franja.CreateFeature(f)

	output.Destroy()

# Calcula el geohash de un punto
def calcular_geohash(lon,lat,nombre,poly_bolivia):
	precision = 9
	hashcode = geohash.encode(lat,lon,9)

	## Generar un SHP con esta franja
	drv = ogr.GetDriverByName("ESRI Shapefile")
	outputfile = "../resultados/geohash_" + nombre + ".shp"
	if os.path.exists(outputfile):
		drv.DeleteDataSource(outputfile)
	output = drv.CreateDataSource(outputfile)
	layer = output.CreateLayer("GeoHash", geom_type=ogr.wkbPolygon)

	## Calcular el espacio con el mismo codigo que el punto
	box = geohash.bbox(hashcode)

	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(box['w'],box['n'])
	ring.AddPoint(box['e'],box['n'])
	ring.AddPoint(box['e'],box['s'])
	ring.AddPoint(box['w'],box['s'])
	ring.AddPoint(box['w'],box['n'])
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	poly_intersection = poly.Clone()
	
	nortesur = ogr.Geometry(ogr.wkbLineString)
	nortesur.AddPoint(box['w'],box['n'])
	nortesur.AddPoint(box['w'],box['s'])
	oesteeste = ogr.Geometry(ogr.wkbLineString)
	oesteeste.AddPoint(box['w'],box['n'])
	oesteeste.AddPoint(box['e'],box['n'])
	
	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
	wgs84 = osr.SpatialReference()
	wgs84.ImportFromEPSG(4326)
	transfLambert = osr.CoordinateTransformation(wgs84,lambert)
	poly_intersection.Transform(transfLambert)
	nortesur.Transform(transfLambert)
	oesteeste.Transform(transfLambert)
	
	print "************************************************"
	print "Codigo geohash - " + nombre + ":"
	print " * valor: " + hashcode
	print " * area equi-codigo:"
	print "   * superficie: %gm2 (equivalente de un cuadrado de %gm de lado)" % (poly_intersection.GetArea(),math.sqrt(poly_intersection.GetArea()))
	print "   * es un rectangulo de perimetro de %gm" % (poly_intersection.Boundary().Length())
	print "                      con tamano norte-sur de %gm" % (nortesur.Length())
	print "                      con tamano oeste-este de %gm" % (oesteeste.Length())
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)

	output.Destroy()

def cargar_capitales():
	# Abrir el archivo
	ds = ogr.Open("../datos/capitales_departamentales.gml")

	# Recuperar la primera y unica capa
	#print ds.GetLayerCount()
	return ds.GetLayer(0)

#DATA_PATH = os.path.dirname("../datos/")
#LIMITES_BOLIVIA_PATH = os.path.join(DATA_PATH, 'limite_nacional.gml')
DATA_PATH = os.path.dirname("../datos/BOL_adm/")
LIMITES_BOLIVIA_PATH = os.path.join(DATA_PATH, 'BOL_adm0.shp')

# Abrir el archivo
ds = ogr.Open(LIMITES_BOLIVIA_PATH)
print "Archivo: %s" % ds.GetName()

# Recuperar la primera y unica capa
layer = ds.GetLayer(0)
print "Capa: %s" % layer.GetName()
origSpatialRef = layer.GetSpatialRef()

# Recuperar la geometria del unico "feature"
feature = layer.GetFeature(0)
poly_bolivia = feature.GetGeometryRef()
#poly_bolivia.CloseRings()
#pol = ogr.Geometry(ogr.wkbPolygon)
#pol.AddGeometry(poly_bolivia)
#print pol.GetGeometryName()

# Definir la proyeccion conica conforme de Lambert de Bolivia
lambert = osr.SpatialReference()
lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Reproyectar en Lambert
SpatialRefTransform = osr.CoordinateTransformation(origSpatialRef,lambert)
poly_bolivia.Transform(SpatialRefTransform)

# Calcular r y tita para cada punto
pts = poly_bolivia.GetGeometryRef(0)
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

print "Valores de x: [%f,%f]" % (min(vx),max(vx))
print "Valores de y: [%f,%f]" % (min(vy),max(vy))
print "Valores de r: [%f,%f]" % (min(vr),max(vr))
print "Valores de tita: [%f,%f]" % (min(vtita),max(vtita))
print "Valores de z: [%d,%d]" % (min(vz),max(vz))

# Numero de codigos
num_codigos = max(vz) - min(vz)
print "Numero de codigos: %d" % num_codigos

# Superficie de Bolivia
area_bolivia_calc = poly_bolivia.GetArea()
area_bolivia_ref = 1098581000000.0
print "Superficie Bolivia (calculado: m2) = %d" % area_bolivia_calc
print "Superficie Bolivia (verdadero: m2) = %d" % area_bolivia_ref

# Superficie por codigo
area_por_codigo = area_bolivia_ref/num_codigos
lado_area_por_codigo = math.sqrt(area_por_codigo)
print "Superficie promedio por codigo (m2) = %d" % area_por_codigo
print "Lado de un cuadrado de la superficie promedio por codigo (m) = %d" % lado_area_por_codigo

# Tests con La Paz, plaza Murillo
#calcular_franja(-68.13352,-16.49568,"La_Paz_plaza_murillo",poly_bolivia)
#calcular_franja(-63.18213,-17.78339,"Santa_Cruz_plaza_24_de_septiembre",poly_bolivia)
#calcular_franja(-66.15696,-17.39377,"Cochabamba_plaza_14_de_septiembre",poly_bolivia)
#calcular_franja(-68.75364,-11.01827,"Cobija_plaza_german_busch",poly_bolivia)
#calcular_franja(-64.73440,-21.53390,"Tarija_plaza_luis_de_fuentes",poly_bolivia)

# Capitales
# Abrir el archivo
ds = ogr.Open("../datos/capitales_departamentales.gml")
layer_capitales = ds.GetLayer(0)

capital = layer_capitales.GetNextFeature()
while capital:
	geom = capital.GetGeometryRef()
	calcular_codigo_ine(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),poly_bolivia)
	calcular_geohash(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),poly_bolivia)
	capital.Destroy()
	capital = layer_capitales.GetNextFeature()
