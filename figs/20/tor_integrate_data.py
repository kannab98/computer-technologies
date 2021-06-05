
import numpy as np
import sys
sys.path.append(".")
from time import time
import pandas as pd

from tasks.task_20 import Tor, rect, monte_carlo

R = (2.5, 0.5)
obj = Tor(a = R)
a = np.sum(R)

N = np.arange(10, 300, 1)
V = np.zeros(N.size)
V0 = np.zeros(N.size)
t = np.zeros(N.size)
t0 = np.zeros(N.size)
V_true = obj.volume

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

d = {"N": N, "V_rect": V, "V_monte": V0, 
     "eps_rect": 100*np.abs(V-V_true)/V_true,
     "eps_monte": 100*np.abs(V0-V_true)/V_true,
     "t_rect": t,
     "t_monte": t0,
}
pd = pd.DataFrame(d)
pd.to_csv("figs/20/tor.csv", sep="\t", index=False)








