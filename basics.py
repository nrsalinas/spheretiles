from math import radians, degrees, sin, cos, tan, acos, asin, atan2
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
