
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.style import use
# use(["science"])

import sys
sys.path.append(".")
from time import time

from tasks.task_20 import Cube, Tor, cart2sphere, rect, monte_carlo, cartesian_meshgrid, cartesian_random, sphere2cart

R = (2.5, 0.5)
obj = Tor(a = R)
a = np.sum(R)



def tor(R=1, r=2):
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


n = 32

r, dr = obj.grid([n,n,n])
r0 = tor(*R)
mask = obj.mask(r)



fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(5,5))

plt.plot(r[0][mask].flat, r[1][mask].flat, r[2][mask].flat, ".")
plt.plot(r0[0].flat, r0[1].flat, r0[2].flat, ".")


plt.show()