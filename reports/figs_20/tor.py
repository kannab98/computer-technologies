import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science"])

def tor(r=1, R=2):
    """
    Генерация тороидальной поверхности
    """
    el = np.linspace(-np.pi, np.pi, 32, endpoint=True)
    az = np.linspace(0, 2*np.pi, 32, endpoint=True)

    el, az = np.meshgrid(az, el)


    x = (R + r*np.cos(el)) * np.cos(az) 
    y = (R + r*np.cos(el)) * np.sin(az) 
    z = r * np.sin(el)

    return x, y, z

import sys
sys.path.append("tasks")

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
