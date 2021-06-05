
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use(["science", "vibrant"])

import pandas as pd

df = pd.read_csv("figs/20/tor.csv", sep="\t", header=0)

N = df.iloc[:,0].values
eps_rect = df.iloc[:,3].values
eps_monte = df.iloc[:,4].values
t_rect = df.iloc[:,5].values
t_monte = df.iloc[:,6].values

plt.figure()
plt.plot(N, eps_rect, label="Monte-Carlo")
plt.plot(N, eps_monte, label="Rectangle")
plt.ylim([None,10])
plt.xlabel("$N$")
plt.ylabel("$\\varepsilon, \\%$")
plt.legend()
plt.savefig("figs/20/eps")

plt.figure()
plt.plot(N, 1e3*t_monte, label="Monte-Carlo")
plt.plot(N, 1e3*t_rect, label="Rectangle")
plt.ylabel("$t, \\text{ms}$")
plt.xlabel("$N$")
plt.legend()
plt.savefig("figs/20/time")



plt.show()




