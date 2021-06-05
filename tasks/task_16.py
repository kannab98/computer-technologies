# Задание 16. Исследуйте движение электрона в неоднородном магнитном поле в зависимости от
# начальных условий и степени однородности поля. Постройте зависимости скорости и
# координаты частицы от времени, а также оцените точность интегрирования в зависимости от
# схемы интегрирования и величины шага интегрирования.
import numpy as np
import matplotlib.pyplot as plt

# Решение систем дифференциальных уравнений
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

q = -1
# q = -1.6021766208*10**(-19)
m = 1
# m_e = 9.1093837015*10**(-31)

B = np.array([1, 0, 0])
V = np.array([0, 0, 1])
r = np.array([0, 1e2, 0])

def func(t, y):
    pos = y[0:3]
    vel = y[3:]
    return np.array([   
                        vel,
                        q/m * np.cross(vel, B)
    ]).flat

t = np.linspace(0, 50, 1024)
sol = solve_ivp(func, t_span=[0, 50], t_eval=t, 
                        y0 = np.hstack((r, V))
)

plt.figure()
plt.plot(sol.t, sol.y[0], label="vx")
plt.plot(sol.t, sol.y[1], label="vy")
plt.plot(sol.t, sol.y[2], label="vz")
plt.legend()
plt.savefig("kek1")

plt.figure()
plt.plot(sol.y[2], sol.y[1])
plt.savefig("kek2")


plt.figure()
plt.plot(sol.t, sol.y[3], label="vx")
plt.plot(sol.t, sol.y[4], label="vy")
plt.plot(sol.t, sol.y[5], label="vz")
plt.legend()
plt.savefig("kek3")