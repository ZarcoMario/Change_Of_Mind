from shapely import Point, MultiPoint
import shapely
import numpy as np
import matplotlib.pyplot as plt

geoms = np.array([Point(0, 0), Point(1, 1), Point(2, 2)])
polygon = shapely.box(0, 0, 2, 2)
res = shapely.contains(polygon, geoms)
x, y = polygon.exterior.xy
plt.plot(x, y)
plt.show()
print(res)