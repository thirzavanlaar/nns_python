import numpy as np
from math import cos,sin,asin,sqrt

# Great circle distance (taking into account curvature of the Earth) and size
def haversine_V(x, y):
   radius = 6371000
   lon1, lat1, size1 = x
   lon2, lat2, size2 = y

   dlat = lat2-lat1
   dlon = lon2-lon1

   a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
   c = 2*asin(sqrt(a))
   distance = radius*c
   V = (size1+size2)/distance
   return V


# Great circle distance (taking into account curvature of the Earth)
def haversine(x, y):
   radius = 6371000
   lon1, lat1 = x
   lon2, lat2 = y

   dlat = lat2-lat1
   dlon = lon2-lon1

   a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
   c = 2*asin(sqrt(a))
   distance = radius*c
   return distance
