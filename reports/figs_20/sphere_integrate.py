
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science"])

import sys
sys.path.append(".")
from time import time

from tasks.task_20 import Sphere, rect, monte_carlo, cartesian_meshgrid, cartesian_random

a = 3
obj = Sphere(a = a)

N = np.arange(10, 200, 5)
V = np.zeros(N.size)
V0 = np.zeros(N.size)
t = np.zeros(N.size)
t0 = np.zeros(N.size)
V_true = obj.volume()

for i, n in enumerate(N):
    r, dr = cartesian_meshgrid([-a, -a, -a], [a, a, a], [n, n, n])
    r0  = cartesian_random([-a, -a, -a], [a, a, a], np.prod([n, n, n]))


    mask = obj.mask(r)
    idx = np.where(mask)

    time1 = time()
    V[i] = rect(mask, dr) 
    time2 = time()
    t[i] = time2- time1
    V0[i] = monte_carlo(mask, r0)
    time3 = time()
    t0[i] = time3- time2

plt.figure()
plt.plot(N, np.abs(V0-V_true)/V_true, label="Monte-Carlo")
plt.plot(N, np.abs(V-V_true)/V_true, label="Rectangle")
plt.legend()

plt.figure()
plt.plot(N, t0, label="Monte-Carlo")
plt.plot(N, t, label="Rectangle")
plt.legend()



plt.show()



