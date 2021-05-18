import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science"])

import sys
sys.path.append("tasks")
from task_20 import tor

x,y,z = tor()

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.plot_wireframe(x,y,z)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
fig.savefig("figs/tor")

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.view_init(azim=0, elev=90)
ax.plot_wireframe(x,y,z)
ax.set_zticks([])
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")


fig.savefig("figs/tor-top")
