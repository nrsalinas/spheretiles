from math import radians, degrees, sin, cos, tan, acos, asin, atan2
from itertools import combinations
import numpy as np

"""
Punta del Este: -54.9330555556, -34.9786111111
Bucaramanga: -73.1161111111, 7.11861111111
Santa Cruz de la Sierra: -63.1975, -17.7891666667
"""


def vectorize(lon, lat):
	x = cos(radians(lon)) * cos(radians(lat))
	y = sin(radians(lon)) * cos(radians(lat))
	z = sin(radians(lat))
	v = np.array([[x,y,z]])
	return v


def coors(vector):
	lon = degrees(atan2(vector[0][1], vector[0][0]))
	lat = degrees(asin(vector[0][2]))
	return (lon, lat)


def dist(point0, point1, radius = 6371.0):
	d = np.dot(point0, point1.transpose())
	return acos(d) * radius


def slerp(point0, point1, fraction):
	w = acos(np.dot(point0, point1.transpose()))
	sl = (sin((1.0-fraction) * w) * point0) / sin(w) \
		+ (sin(fraction * w) * point1) / sin(w)
	return sl

def icosahedron_vertices():
	h = 1.0 / (5.0)**0.5
	xminus = 0.5 - ((5.0)**0.5 / 10.0)
	xplus = 0.5 + ((5.0)**0.5 / 10.0)
	vertices = [
		np.array([[0.0, 0.0, 1.0]]),
		np.array([[0.0, 0.0, -1.0]]),

		np.array([[-2*h, 0.0, h]]),
		np.array([[2*h, 0.0, -h]]),

		np.array([[ xplus, (xminus ** 0.5), h]]),
		np.array([[ -xplus, -(xminus ** 0.5), -h]]),

		np.array([[ -xminus, -(xplus ** 0.5), h]]),
		np.array([[ xminus, (xplus ** 0.5), -h]]),

		np.array([[ -xminus, (xplus ** 0.5), h]]),
		np.array([[ xminus, -(xplus ** 0.5), -h]]),

		np.array([[ xplus, -(xminus ** 0.5), h]]),
		np.array([[ -xplus, (xminus ** 0.5), -h]])
		]

	return vertices

class polyhedron:

	def __init__(self):
		self.vertices = icosahedron_vertices()  # list of Numpy arrays
		self.neighboors = {0: [2, 4, 6, 8, 10],
						 1: [3, 5, 7, 9, 11],
						 2: [0, 5, 6, 8, 11],
						 3: [1, 4, 7, 9, 10],
						 4: [0, 3, 7, 8, 10],
						 5: [1, 2, 6, 9, 11],
						 6: [0, 2, 5, 9, 10],
						 7: [1, 3, 4, 8, 11],
						 8: [0, 2, 4, 7, 11],
						 9: [1, 3, 5, 6, 10],
						 10: [0, 3, 4, 6, 9],
						 11: [1, 2, 5, 7, 8]}
		self.triangles = []
		self.edges = []
		for v0 in self.neighboors:
			for v1 in self.neighboors[v0]:
				if v0 < v1:
					self.edges.append((v0,v1))


	def set_triangles(self):
		self.triangles = []
		for vi in self.neighboors:
			for ni0,ni1 in combinations(self.neighboors[vi], 2):
				if ni0 in self.neighboors[ni1] and ni1 in self.neighboors[ni0]:
					this_triangle = sorted([vi, ni0, ni1])
					if this_triangle not in self.triangles:
						self.triangles.append(this_triangle)
		return None

	def bisect_edges(self):
		for vi0,vi1 in self.edges:
			new_vertex = slerp(self.vertices[vi0],self.vertices[vi1],0.5)
			newvi = len(self.vertices)
			self.vertices.append(new_vertex)
			self.neighboors[newvi] = [vi0, vi1]
			self.neighboors[vi0].append(newvi)
			self.neighboors[vi0].remove(vi1)
			self.neighboors[vi1].append(newvi)
			self.neighboors[vi1].remove(vi0)

		# get remaining neoghboors of new vertices
		for vi in self.neighboors:
			if len(self.neighboors[vi]) == 2:
				neigh0 = self.neighboors[vi][0]
				neigh1 = self.neighboors[vi][1]
				for conneigh0, conneigh1 in [(x,y) for x in self.neighboors[neigh0]
							for y in self.neighboors[neigh1] if x != vi and y != vi]:
					if 0 < dist(self.vertices[conneigh0], self.vertices[conneigh1]) \
							< dist(self.vertices[vi], self.vertices[neigh0]) * 1.2 \
							and conneigh0 not in self.neighboors[vi] and conneigh1 \
							not in self.neighboors[vi]:
						self.neighboors[vi].append(conneigh0)
						self.neighboors[vi].append(conneigh1)
					if len(self.neighboors[vi]) == 6:
						self.neighboors[vi].sort()
						break

		self.edges = []
		for v0 in self.neighboors:
			for v1 in self.neighboors[v0]:
				if v0 < v1:
					self.edges.append((v0,v1))

		return None

	def find_closest_vertex(self, vector):
		advance = True
		tvertex = 0
		tdist = dist(self.vertices[tvertex], vector)
		while advance:
			advance = False
			for tneigh in self.neighboors[tvertex]:
				dist2neigh = dist(self.vertices[tneigh], vector)
				if dist2neigh < tdist:
					tdist = dist2neigh
					tvertex = tneigh
					advance = True
		return tvertex


	def vert_markers(self):
		out = ""
		for indx, vector in enumerate(self.vertices):
			lon, lat = coors(vector)
			out += "\t\tvar marker{0} = WE.marker([{1}, {2}]).addTo(earth);\n".format(indx, lat, lon)
			out += "\t\tmarker%s.bindPopup(\"<b>%s : %.2f, %.2f</b><br>%s</br>\", {maxWidth: 150})\n" % (indx, indx, lat, lon, self.neighboors[indx])
		return out

	def polygonize(self):
		self.set_triangles()
		out = "\n"
		for it, tri in enumerate(self.triangles):
			out += "\t\tvar triangle_{0} = WE.polygon([".format(it)
			myvectors = map(lambda x: self.vertices[x], tri)
			#print 'myvectors', myvectors
			mycoors = map(lambda x: coors(x), myvectors)
			#print 'mycoors', mycoors
			out += ",".join(map(lambda x: "[%s , %s]"  % (x[1], x[0]), mycoors))
			out += "]).addTo(earth);\n"

		return out

	def to_html(self, polygons = True, vertices = False):
		bffr = """
<!DOCTYPE HTML>
<html>
  <head>
    <script src="http://www.webglearth.com/v2/api.js"></script>
    <script>
      function initialize() {
        var earth = new WE.map('earth_div');

        WE.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',{
          attribution: 'Copyright OpenStreetMap contributors'
        }).addTo(earth);"""

		if vertices:
			bffr += self.vert_markers() + "\n"

		if polygons:
			bffr += self.polygonize()

		bffr += """
		}
    </script>
    <style>
      html, body{padding: 0; margin: 0;}
      #earth_div{top: 0; right: 0; bottom: 0; left: 0; position: absolute !important;}
    </style>
    <title>WebGL Earth API: Geodesic grid.</title>
  </head>
  <body onload="initialize()">
    <div id="earth_div"></div>
  </body>
</html>
		"""
		with open('index.html','w') as fhandle:
			fhandle.write(bffr)
