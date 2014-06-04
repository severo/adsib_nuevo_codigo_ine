# Para hacer funcionar con Debian
#  sudo aptitude install python-geohash python-gdal python-pip
# instalar las librerias python mgrs y mlocs:
#  sudo pip install mgrs
#  sudo pip install mlocs
# Hay un problema con la instalacion de mlocs, para terminar:
#  sudo mv build/mlocs/Library/Python/2.7/site-packages/mlocs /usr/local/lib/python2.7/dist-packages/
#  sudo chown root:staff /usr/local/lib/python2.7/dist-packages/mlocs
#  sudo rm -rf build

import csv
import geohash
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import mgrs
import mlocs
import numpy as np
from osgeo import ogr
import os
import osr

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
def codigo_ine(lon,lat,nombre,num_carac,simple=False):

	if num_carac != 9:
		if simple:
			return ''
		else:
			return '',float('NaN'),float('NaN'),float('NaN')

	(r,tita,z) = codigo_ine_wgs84(lon,lat)

	if simple:
		return r


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

	poly_intersection = poly.Intersection(poly_bolivia_lambert)

	superficie = poly_intersection.GetArea()
	largo = poly_intersection.Boundary().Length()/2
	ancho = float(sum(vdeltar)) / len(vdeltar)
	ancho_min = min(vdeltar)
	ancho_max = max(vdeltar)
	print "************************************************"
	print "Codigo INE z - " + nombre + ":"
	print " * valor: %d" % z
	print " * area equi-codigo: arco de %gm2, %gm de largo por %gm de ancho promedio (min=%gm max=%gm)" % (superficie,largo,ancho,ancho_min,ancho_max)
	print "************************************************"

	poly_intersection.Transform(transfWGS84)
	f = ogr.Feature(feature_def=franja.GetLayerDefn())
	f.SetGeometryDirectly(poly_intersection)
	franja.CreateFeature(f)
	output.Destroy()

	return z,superficie,largo,ancho

# Calcula el geohash de un punto
def codigo_geohash(lon,lat,nombre,num_carac,simple=False):
	hashcode = geohash.encode(lat,lon,num_carac)

	if simple:
		return hashcode

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

	superficie = poly_intersection.GetArea()
	ns = nortesur.Length()
	oe = oesteeste.Length()
	print "************************************************"
	print "Codigo geohash - " + nombre + ":"
	print " * valor: " + hashcode
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (superficie,ns,oe)
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)
	output.Destroy()

	return hashcode,superficie,ns,oe

# Calcula el codigo MGRS de un punto
def codigo_mgrs(lon,lat,nombre,num_carac,simple=False):
	m = mgrs.MGRS()

	tabla_precision = {5:0,7:1,9:2,11:3,13:4,15:5}
	# Ver http://en.wikipedia.org/wiki/Military_grid_reference_system
	tabla_km = {5:100000,7:10000,9:1000,11:100,13:10,15:1}

	if num_carac not in tabla_precision:
		if simple:
			return ''
		else:
			return '',float('NaN'),float('NaN'),float('NaN')

	# Nota: ponemos precision=2 para llegar a un codigo de 9 caracteres
	mgrscode = m.toMGRS(lat,lon,True,tabla_precision[num_carac])

	if simple:
		return mgrscode

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
	prec = tabla_km[num_carac]/2

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

	superficie = poly_intersection.GetArea()
	ns = nortesur.Length()
	oe = oesteeste.Length()
	print "************************************************"
	print "Codigo mgrs - " + nombre + ":"
	print " * valor: " + mgrscode
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (superficie,ns,oe)
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	poly.Transform(transfWGS84)
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)
	output.Destroy()

	return mgrscode,superficie,ns,oe
	
# Calcula el codigo MLOCS (Maidenhead Locator System) de un punto
def codigo_mlocs(lon,lat,nombre,num_carac,simple=False):

	tabla_precision = {2:1,4:2,6:3,8:4,10:5,12:6,14:7}
	# Ver http://en.wikipedia.org/wiki/Maidenhead_Locator_System#Description_of_the_system
	tabla_grados_lat = {2:10.0, 4:10.0/10, 6:10.0/(10*24), 8:10.0/(10*24*10), 10:10.0/(10*24*10*24), 12:10.0/(10*24*10*24*10), 14:10.0/(10*24*10*24*10*24)}
	tabla_grados_lon = {2:20.0, 4:20.0/10, 6:20.0/(10*24), 8:20.0/(10*24*10), 10:20.0/(10*24*10*24), 12:20.0/(10*24*10*24*10), 14:10.0/(10*24*10*24*10*24)}

	if num_carac not in tabla_precision:
		if simple:
			return ''
		else:
			return '',float('NaN'),float('NaN'),float('NaN')

	mlocscode = mlocs.toMaiden([lat,lon],tabla_precision[num_carac])

	if simple:
		return mlocscode

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
	precx = tabla_grados_lon[num_carac]
	precy = tabla_grados_lat[num_carac]

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

	superficie = poly_intersection.GetArea()
	ns = nortesur.Length()
	oe = oesteeste.Length()
	print "************************************************"
	print "Codigo mlocs - " + nombre + ":"
	print " * valor: " + mlocscode
	print " * area equi-codigo: rectangulo de %gm2, lado NS: %gm, OE: %gm" % (superficie,ns,oe)
	print "************************************************"

	f = ogr.Feature(feature_def=layer.GetLayerDefn())
	f.SetGeometryDirectly(poly)
	layer.CreateFeature(f)
	output.Destroy()

	return mlocscode,superficie,ns,oe

def codigo(algo,lon,lat,nombre,num_carac,simple=False):
	if algo == "geohash":
		return codigo_geohash(lon,lat,nombre,num_carac,simple)
	elif algo == "mgrs":
		return codigo_mgrs(lon,lat,nombre,num_carac,simple)
	elif algo == "mlocs":
		return codigo_mlocs(lon,lat,nombre,num_carac,simple)
	elif algo == "ine":
		return codigo_ine(lon,lat,nombre,num_carac,simple)

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
	codigo_ine(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),9)
	codigo_geohash(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),9)
	codigo_mgrs(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),9)
	codigo_mlocs(geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),8)

	if capital.GetFieldAsString("NOMBRE") == "LA PAZ":
		code = {'geohash': [],'mgrs': [],'mlocs': [],'ine': []}
		superficie = {'geohash': [],'mgrs': [],'mlocs': [],'ine': []}
		precision = {'geohash': [],'mgrs': [],'mlocs': [],'ine': []}
		nortesur = {'geohash': [],'mgrs': [],'mlocs': [],'ine': []}
		oesteeste = {'geohash': [],'mgrs': [],'mlocs': [],'ine': []}
		for algo in ['geohash','mgrs','mlocs','ine']:
			for num_carac in range(4,14):
				c,s,n,o = codigo(algo,geom.GetX(),geom.GetY(),capital.GetFieldAsString("NOMBRE"),num_carac)
				code[algo].append(c)
				superficie[algo].append(s)
				precision[algo].append(math.sqrt(s))
				nortesur[algo].append(n)
				oesteeste[algo].append(o)
	capital.Destroy()
	capital = layer_capitales.GetNextFeature()

## Graficos

def escala_sistema_metrico(l,pos):
	if l >= 1000:
		return "%d km" % (l/1000)
	elif l >= 1:
		return "%d m" % l
	elif l >= 0.001:
		return "%d mm" % (l*1000)
	else:
		return "%f mm" % l

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16,
        }


x = range(4,14)
y_geohash = precision['geohash']
y_mgrs = precision['mgrs']
y_mlocs = precision['mlocs']
y_ine = precision['ine']

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
#ax.yaxis.set_major_formatter(tick.ScalarFormatter())
#ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(tick.FuncFormatter(escala_sistema_metrico))

plt.plot(x, y_ine, 'kD',x, y_geohash, 'ko',x, y_mgrs, 'k^',x, y_mlocs, 'k*')
plt.title('Precision del codigo segun caracteres (punto en La Paz)', fontdict=font)
#plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)
plt.xlabel('Numero de caracteres', fontdict=font)
plt.ylabel('Precision', fontdict=font)
#plt.yscale('log')
# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.xlim(3.5,13.5)
plt.legend((r'INE',r'geohash', r'mgrs', r'mlocs'))

plt.grid(True)

#plt.show()
plt.savefig("../resultados/precision_codigos.png",format="png")

## Busqueda de duplicados

# cargar los dos archivos csv de manzanas y localidades

code = {}
dup = {}

for algo in ['ine','geohash','mgrs','mlocs']:
	print algo
	if algo not in code:
		code[algo] = {}
	if algo not in dup:
		dup[algo] = {}

	for num_carac in range(0,14):
		print num_carac
		for row in csv.DictReader(file('../datos/INE/localidades.csv'),delimiter=','):
			lat = float(row['Y'])
			lon = float(row['X'])

			c = codigo(algo,lon,lat,"test",num_carac,True)

			if num_carac not in code[algo]:
				code[algo][num_carac] = {}
			if num_carac not in dup[algo]:
				dup[algo][num_carac] = 0

			if c:
				if c in code[algo][num_carac]:
					dup[algo][num_carac] += 1
				code[algo][num_carac][c] = code[algo][num_carac].get(c,0) + 1

		for row in csv.DictReader(file('../datos/INE/manzanas.csv'),delimiter=','):
			lat = float(row['lat'])
			lon = float(row['long'])

			c = codigo(algo,lon,lat,"test",num_carac,True)

			if num_carac not in code[algo]:
				code[algo][num_carac] = {}
			if num_carac not in dup[algo]:
				dup[algo][num_carac] = 0

			if c:
				if c in code[algo][num_carac]:
					dup[algo][num_carac] += 1
				code[algo][num_carac][c] = code[algo][num_carac].get(c,0) + 1

		# buscar duplicados
		print "len %d" % len(code[algo][num_carac])
		print "duplicados: %d" % dup[algo][num_carac]
