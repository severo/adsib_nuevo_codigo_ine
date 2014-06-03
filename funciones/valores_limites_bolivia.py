# Para hacer funcionar con Debian
#  sudo aptitude install python-geohash python-gdal python-pip
# instalar las librerias python mgrs y mlocs:
#  sudo pip install mgrs
#  sudo pip install mlocs

import os
from osgeo import ogr
import osr
import geohash
import math
import mgrs
import mlocs
import numpy as np

def transf_wgs84_lambert(sentido=0):
	# WGS84
	wgs84 = osr.SpatialReference()
	wgs84.ImportFromEPSG(4326)
	# Proyeccion conica conforme de Lambert de Bolivia
	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
	# Tranformaciones
	if (sentido==0):
		return osr.CoordinateTransformation(wgs84,lambert)
	else:
		return osr.CoordinateTransformation(lambert,wgs84)

def wgs84_to_lambert():
	return transf_wgs84_lambert(sentido=0)

def lambert_to_wgs84():
	return transf_wgs84_lambert(sentido=1)
	
def codigo_ine_lambert(x,y):
	# Calcular r, tita y z
	r = math.sqrt(x**2+ y**2)
	tita = np.arctan(y/x)+np.pi
	z = math.floor(r * math.exp(tita))

	return (r,tita,z)

def codigo_ine_wgs84(lon,lat):
	x,y,z = wgs84_to_lambert().TransformPoint(lon, lat, 0)
	
	return codigo_ine_lambert(x,y)

def codigos_ine_bolivia(poly_bolivia_lambert):
	# Calcular r y tita para cada punto
	pts = poly_bolivia_lambert.GetGeometryRef(0)
	vx = []
	vy = []
	vr = []
	vtita = []
	vz = []
	for p in xrange(pts.GetPointCount()):
		x = pts.GetX(p)
		y = pts.GetY(p)
		r,tita,z = codigo_ine_lambert(x,y)
		vx.append(x)
		vy.append(y)
		vr.append(r)
		vtita.append(tita)
		vz.append(z)

	# Min y Max
	print "************************************************"
	print "Analisis del codigo z para Bolivia"
	print " * Valores de x: [%g,%g]" % (min(vx),max(vx))
	print " * Valores de y: [%g,%g]" % (min(vy),max(vy))
	print " * Valores de r: [%g,%g]" % (min(vr),max(vr))
	print " * Valores de tita: [%g,%g]" % (min(vtita),max(vtita))
	print " * Valores de z: [%d,%d]" % (min(vz),max(vz))

	# Numero de codigos
	num_codigos = max(vz) - min(vz)
	print " * Numero de codigos: %d" % num_codigos

	# Superficie de Bolivia
	area_bolivia_calc = poly_bolivia_lambert.GetArea()
	area_bolivia_ref = 1098581000000.0
	print " * Superficie Bolivia: calculado: %gm2, oficial: %gm2" % (area_bolivia_calc,area_bolivia_ref)

	# Superficie por codigo
	area_por_codigo = area_bolivia_ref/num_codigos
	lado_area_por_codigo = math.sqrt(area_por_codigo)
	print " * Superficie promedio por codigo (m2) = %g" % area_por_codigo
	print " * Lado de un cuadrado de la superficie promedio por codigo (m) = %g" % lado_area_por_codigo
	
	return (vx,vy,vr,vtita,vz)

# Calcula el codigo INE que pasa por un punto
def calcular_codigo_ine(lon,lat,nombre,poly_bolivia,vtita):
	(r,tita,z) = codigo_ine_wgs84(lon,lat)

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
	print " * area equi-codigo: arco de %gm2, %gm de largo por %gm de ancho promedio (min=%gm max=%gm)" % (poly_intersection.GetArea(),poly_intersection.Boundary().Length()/2,float( sum(vdeltar) ) / len(vdeltar),min(vdeltar),max(vdeltar))
	print "************************************************"
	
	poly_intersection.Transform(transfWGS84)
	f = ogr.Feature(feature_def=franja.GetLayerDefn())
	f.SetGeometryDirectly(poly_intersection)
	franja.CreateFeature(f)

	output.Destroy()

# Calcula el geohash de un punto
def calcular_geohash(lon,lat,nombre):
	precision = 9
	hashcode = geohash.encode(lat,lon,9)

	## Generar un SHP
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
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (poly_intersection.GetArea(),nortesur.Length(),oesteeste.Length())
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)

	output.Destroy()

# Calcula el codigo MGRS de un punto
def calcular_mgrs(lon,lat,nombre):
	m = mgrs.MGRS()
	# Nota: ponemos precision=2 para llegar a un codigo de 9 caracteres
	mgrscode = m.toMGRS(lat,lon,True,2)
		
	centro = m.toLatLon(mgrscode)
	pt_centro = ogr.Geometry(ogr.wkbPoint)
	pt_centro.AddPoint(centro[1], centro[0])
	
	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
	wgs84 = osr.SpatialReference()
	wgs84.ImportFromEPSG(4326)
	transfLambert = osr.CoordinateTransformation(wgs84,lambert)
	transfWGS84 = osr.CoordinateTransformation(lambert,wgs84)

	pt_centro.Transform(transfLambert)
	x_centro = pt_centro.GetX()
	y_centro = pt_centro.GetY()
	
	# precision = 2 significa precision de 1km
	prec = 500
	
	## Calcular el espacio con el mismo codigo que el punto
	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(x_centro-prec,y_centro-prec)
	ring.AddPoint(x_centro+prec,y_centro-prec)
	ring.AddPoint(x_centro+prec,y_centro+prec)
	ring.AddPoint(x_centro-prec,y_centro+prec)
	ring.AddPoint(x_centro-prec,y_centro-prec)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	poly_intersection = poly.Clone()
	
	nortesur = ogr.Geometry(ogr.wkbLineString)
	nortesur.AddPoint(x_centro+prec,y_centro-prec)
	nortesur.AddPoint(x_centro+prec,y_centro+prec)
	oesteeste = ogr.Geometry(ogr.wkbLineString)
	oesteeste.AddPoint(x_centro-prec,y_centro-prec)
	oesteeste.AddPoint(x_centro+prec,y_centro-prec)
	
	## Generar un SHP
	drv = ogr.GetDriverByName("ESRI Shapefile")
	outputfile = "../resultados/mgrs_" + nombre + ".shp"
	if os.path.exists(outputfile):
		drv.DeleteDataSource(outputfile)
	output = drv.CreateDataSource(outputfile)
	layer = output.CreateLayer("MGRS", geom_type=ogr.wkbPolygon)

	print "************************************************"
	print "Codigo mgrs - " + nombre + ":"
	print " * valor: " + mgrscode
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (poly_intersection.GetArea(),nortesur.Length(),oesteeste.Length())
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	poly.Transform(transfWGS84)
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)
	
	output.Destroy()

# Calcula el codigo MLOCS (Maidenhead Locator System) de un punto
def calcular_mlocs(lon,lat,nombre):
	# Nota: ponemos precision=4 para llegar a un codigo de 8 caracteres
	mlocscode = mlocs.toMaiden([lat,lon],4)

	centro = mlocs.toLoc(mlocscode)
	x_centro = centro[1]
	y_centro = centro[0]
	
	lambert = osr.SpatialReference()
	lambert.ImportFromProj4("+proj=lcc +lat_1=-11.5 +lat_2=-21.5 +lat_0=-24 +lon_0=-64 +x_0=1000000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
	wgs84 = osr.SpatialReference()
	wgs84.ImportFromEPSG(4326)
	transfLambert = osr.CoordinateTransformation(wgs84,lambert)
	transfWGS84 = osr.CoordinateTransformation(lambert,wgs84)
	
	# precision = 4 significa precision de 15" en lat y 30" en lon
	precx = 1.0/240
	precy = 1.0/120
	
	## Calcular el espacio con el mismo codigo que el punto
	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(x_centro-precx,y_centro-precy)
	ring.AddPoint(x_centro+precx,y_centro-precy)
	ring.AddPoint(x_centro+precx,y_centro+precy)
	ring.AddPoint(x_centro-precx,y_centro+precy)
	ring.AddPoint(x_centro-precx,y_centro-precy)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	poly_intersection = poly.Clone()
	
	nortesur = ogr.Geometry(ogr.wkbLineString)
	nortesur.AddPoint(x_centro+precx,y_centro-precy)
	nortesur.AddPoint(x_centro+precx,y_centro+precy)
	oesteeste = ogr.Geometry(ogr.wkbLineString)
	oesteeste.AddPoint(x_centro-precx,y_centro-precy)
	oesteeste.AddPoint(x_centro+precx,y_centro-precy)

	poly_intersection.Transform(transfLambert)
	nortesur.Transform(transfLambert)
	oesteeste.Transform(transfLambert)
	
	## Generar un SHP
	drv = ogr.GetDriverByName("ESRI Shapefile")
	outputfile = "../resultados/mlocs_" + nombre + ".shp"
	if os.path.exists(outputfile):
		drv.DeleteDataSource(outputfile)
	output = drv.CreateDataSource(outputfile)
	layer = output.CreateLayer("MLOCS", geom_type=ogr.wkbPolygon)

	print "************************************************"
	print "Codigo mlocs - " + nombre + ":"
	print " * valor: " + mlocscode
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (poly_intersection.GetArea(),nortesur.Length(),oesteeste.Length())
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)
	
	output.Destroy()

def cargar_capitales():
	# Abrir el archivo
	ds = ogr.Open("../datos/capitales_departamentales.gml")

	# Recuperar la primera y unica capa
	return ds.GetLayer(0)

def cargar_poly_bolivia_lambert():
	DATA_PATH = os.path.dirname("../datos/BOL_adm/")
	LIMITES_BOLIVIA_PATH = os.path.join(DATA_PATH, 'BOL_adm0.shp')

	# Abrir el archivo
	ds = ogr.Open(LIMITES_BOLIVIA_PATH)

	# Recuperar la primera y unica capa
	layer = ds.GetLayer(0)

	# Recuperar la geometria del unico "feature"
	feature = layer.GetFeature(0)
	poly_bolivia = feature.GetGeometryRef()

	# Reproyectar en Lambert
	poly_bolivia.Transform(wgs84_to_lambert())

	return poly_bolivia.Clone()


poly_bolivia_lambert = cargar_poly_bolivia_lambert()

vx,vy,vr,vtita,vz = codigos_ine_bolivia(poly_bolivia_lambert)

# Capitales
# Abrir el archivo
ds = ogr.Open("../datos/capitales_departamentales.gml")
layer_capitales = ds.GetLayer(0)

capital = layer_capitales.GetNextFeature()
while capital:
	geom = capital.GetGeometryRef()
	calcular_codigo_ine(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),poly_bolivia_lambert,vtita)
	calcular_geohash(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"))
	calcular_mgrs(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"))
	calcular_mlocs(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"))
	capital.Destroy()
	capital = layer_capitales.GetNextFeature()
