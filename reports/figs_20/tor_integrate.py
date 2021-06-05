
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science"])

import sys
sys.path.append(".")
from time import time

from tasks.task_20 import Tor, rect, monte_carlo, cartesian_meshgrid, cartesian_random

R = (2.5, 0.5)
obj = Tor(a = R)
a = np.sum(R)

N = np.arange(10, 200, 5)
V = np.zeros(N.size)
V0 = np.zeros(N.size)
t = np.zeros(N.size)
t0 = np.zeros(N.size)
V_true = obj.volume()

for i, n in enumerate(N):
    r, dr = obj.grid([n,n,n], dtype="linear")
    r0  = obj.grid([n,n,n], dtype="random")


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
plt.plot(N, 100*np.abs(V0-V_true)/V_true, label="Monte-Carlo")
plt.plot(N, 100*np.abs(V-V_true)/V_true, label="Rectangle")
plt.xlabel("$N$")
plt.ylabel("$\\varepsilon, \\%$")
plt.legend()

plt.figure()
plt.plot(N, 1e3*t0, label="Monte-Carlo")
plt.plot(N, 1e3*t, label="Rectangle")
plt.ylabel("$t, \\text{ms}$")
plt.xlabel("$N$")
plt.legend()



plt.show()




