import numpy as np
from math import cos,sin,asin,sqrt

# Great circle distance (taking into account curvature of the Earth) and size
def euclidean_V(x, y):
    radius = 6371000
    x1, y1, size1 = x
    x2, y2, size2 = y
 
    part1 = (x2-x1)**2
    part2 = (y2-y1)**2
    distance = sqrt(part1+part2)

    V = (size1+size2)/distance
    return V


def euclidean(x, y):
    radius = 6371000
    x1, y1 = x
    x2, y2 = y
 
    part1 = (x2-x1)**2
    part2 = (y2-y1)**2
    distance = sqrt(part1+part2)

    return distance
