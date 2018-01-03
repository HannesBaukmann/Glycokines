import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

xyz = []
cur = [0, 0, 0]

for _ in range(20):
    axis = random.randrange(0, 3)
    cur[axis] += random.choice([-1, 1])
    xyz.append(cur[:])

x, y, z = zip(*xyz)
ax.plot(x, y, z, label='Random walk')
ax.scatter(x[-1], y[-1], z[-1], c='b', marker='o')   # End point
ax.legend()
plt.show()