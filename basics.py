from math import radians, degrees, sin, cos, tan, acos, asin, atan2
import numpy as np

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
