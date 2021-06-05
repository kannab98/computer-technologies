
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science"])

import sys
sys.path.append(".")

from tasks.task_20 import obj, rect, monte_carlo, cartesian_meshgrid

a = 2
r, dr = cartesian_meshgrid([-a, -a, -a], [a, a, a], [50, 50, 50])
obj = obj(a = a)
mask = obj.mask(r)
idx = np.where(mask)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(r[0][idx], r[1][idx], r[2][idx])

ax.set_xlim(-a,a)
ax.set_ylim(-a,a)
ax.set_zlim(-a,a)
plt.savefig("reports/figs_20/obj_scatter.png", dpi=300)