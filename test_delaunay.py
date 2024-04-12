import numpy as np
points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
from scipy.spatial import Delaunay
tri = Delaunay(points)
p = np.array([[0.1, 0.2], [1.5, 0.5], [0.5, 1.05]])
res = tri.find_simplex(p) < 0
print(res)
import matplotlib.pyplot as plt
plt.triplot(points[:,0], points[:,1], tri.simplices)
plt.plot(points[:,0], points[:,1], 'o')
for p_ in p:
    print(p_)
    plt.plot(p_[0],p_[1], 'ro')
# plt.plot(0.8, 0.2, 'ro')
plt.show()