import numpy as np

# Решение систем уравнений
from scipy.optimize import fsolve
# Решение систем дифференциальных уравнений
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def predator_prey(t: np.array, z: np.ndarray, p:np.ndarray, dtype="LV"):
    x, y = z
    if dtype == "LV":
        func = [
            p[0][0]*x - p[0][1]*x*y,
            p[1][0]*x*y - p[1][1]*y
        ]
    elif dtype == "McA":
        func = [
            p[0][0]*x - p[0][1]*x*y - p[0][2]*x*x,
            p[1][0]*x*y - p[1][1]*y - p[1][2]*y*y,
        ]
    return func

alpha = 1.1
beta = 0.4
gamma = 1.1
delta = 0.4



p = np.array([[alpha, beta, 0], [gamma, delta, 1]])


# Changed integrator. This one gives a more correct solution
# t_span -- Interval of integration (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
# y0 -- Initial state. [prey, predator]

sol = solve_ivp(predator_prey, t_span=[0, 500], y0=[10, 2], args=(p,"McA"))


plt.figure()
plt.plot(sol.y[0])
plt.ylabel("Зайцы")
plt.xlabel("Время")


plt.figure()
plt.plot(sol.y[1])
plt.ylabel("Волки")
plt.xlabel("Время")

plt.figure()
plt.plot(sol.y[0], sol.y[1], ".-")
plt.xlabel("Зайцы")
plt.ylabel("Волки")


plt.show()
