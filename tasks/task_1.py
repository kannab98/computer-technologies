import numpy as np

# Решение систем уравнений
from scipy.optimize import fsolve
# Решение систем дифференциальных уравнений
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt



def dx(t, x, y, a, b):
    return (a - b*x) * y

def dy(t, y, x, a, b):
    return (-a + b*x)

def predator_prey(x, a, b, c, d):
    t = [0, 10]
    return [
        solve_ivp(dx, t_span=t, y0=[x[0]], args=(x[1], a, b)).y,
        solve_ivp(dy, t_span=t, y0=[x[1]], args=(x[0], c, d)).y
    ]

root = predator_prey([50, 100], 1,1,1,1)

plt.plot(root[0][0])
plt.savefig("kek.png")



