import numpy as np
# Решение систем дифференциальных уравнений
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from matplotlib.style import use
# use(["ieee"])

# Разработайте алгоритм решения задачи для модели, описывающей вынужденные
# колебания системы, включающей три тела, соединенные пружинами, при наличии силы вязкого
# сопротивления. Постройте осциллограммы и фазовые портреты колебаний, а также оцените
# точность интегрирования в зависимости от схемы интегрирования и величины шага
# интегрирования.

def func(t, y, d, Q, F):
    return np.array([
        y[3],
        y[4],
        y[5],
        - (1 + d[0]) * y[0] + (1+d[0])/2 * y[1] - 1/Q[0] * y[3],
        + (1 + d[1])/2 * y[0] - (1+d[1]) * y[1] + (1+d[1])/2 * y[2] - 1/Q[1] * y[4] + F*np.cos(t),
        + (1 + d[2])/2 * y[1] - (1+d[2]) * y[2] - 1/Q[2] * y[5]
    ])


y0 = [0, 0, 0, 0, 0, 0]
d = [1, 0, 2]
Q = [1, 1000, 1]
F = 1

t = np.linspace(0, 4*np.pi, 1024)
sol = solve_ivp(func, t_span=[t.min(), t.max()], t_eval=t, 
                        y0 = y0, args=(d, Q, F) )



idx = [1, 2, 3]

# Осциллограммы скорости
plt.figure()
N = [3, 4, 5]
for i, n in enumerate(N):
    plt.plot(t, sol.y[n], "-", label="$\\dot x_{%d}$" % idx[i])
plt.xlabel("$\\tau, \\text{rad}$")
plt.legend()

# Осциллограммы координаты
plt.figure()
N = [0, 1, 2]
for i, n in enumerate(N):
    plt.plot(t, sol.y[n], "-", label="$x_{%d}$" % idx[i])
plt.xlabel("$\\tau, \\text{rad}$")
plt.legend()


# Фазовые портреты
plt.figure()
N = [0, 1, 2]
for i, n in enumerate(N):
    plt.plot(sol.y[n], sol.y[n+3], "-", label="$x_{%d}$" % idx[i])
plt.xlabel("$x$")
plt.ylabel("$\\dot x$")
plt.legend()

plt.show()