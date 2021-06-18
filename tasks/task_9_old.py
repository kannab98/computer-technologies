# Задание 9. Рассчитайте эффективную траекторию ракеты, предназначенной для наименее 
# затратного по топливу запуска с Земли искусственного спутника Марса. Постройте зависимости 
# скорости и координаты ракеты от времени, а также оцените точность интегрирования в 
# зависимости от схемы интегрирования и величины шага интегрирования

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

G = 6.67430*10**(-11)
M_SUN = 1.989*10**30 # kg

M_MARS = 6.39*10**23 # kg
R_MARS = 2.06655*10**11 # m
V_MARS = 24130 # m/s
ANGLE_V_MARS = V_MARS/R_MARS # rad/s
OFFSET_MARS = 0.29*np.pi

M_EARTH = 5.972*10**24 # kg
R_EARTH = 150*10**9 # m
V_EARTH = 29783 # m/s
ANGLE_V_EARTH = V_EARTH/R_EARTH # rad/s

M_ROCKET = 3*10**6

t1 = 0
t2 = 1.9*10**8
t_dV = [t1, t2]
t_acc = 1 # s
Vp = np.sqrt((V_MARS**2 + V_EARTH**2)/2)
dV1 = V_EARTH*(V_EARTH/Vp - 1)
dV2 = V_MARS*(1 - V_MARS/Vp)
dV = [dV1, dV2]

def acc_dV(t, num):
    acc = dV[num] * 1/t_acc * (np.heaviside(t-t_dV[num], 1) - np.heaviside(t-t_dV[num]-t_acc, 1))
    if num == 0:
        acc = np.array([0, acc])
    else:
        acc = np.array([0, -acc])
    return acc

def acc_g(t, r):
    x = r[0]
    y = r[1]
    x_E = R_EARTH*np.cos(ANGLE_V_EARTH*t)
    y_E = R_EARTH*np.sin(ANGLE_V_EARTH*t)
    x_M = R_MARS*np.cos(ANGLE_V_MARS*t + OFFSET_MARS)
    y_M = R_MARS*np.sin(ANGLE_V_MARS*t + OFFSET_MARS)

    R2E = np.sqrt((x-x_E)**2 + (y-y_E)**2)
    R2M = np.sqrt((x-x_M)**2 + (y-y_M)**2)
    R2S = np.sqrt(x**2 + y**2)

    acc_E = M_EARTH/R2E**3 * np.array([x_E - x, y_E - y])
    acc_M = M_MARS/R2M**3 * np.array([x_M - x, y_M - y])
    acc_S = -M_SUN/R2S**3 * np.array([x, y])
    
    # acc = G*(acc_S)
    # acc = G*(acc_M + acc_S)
    acc = G*(acc_E + acc_M + acc_S)
    return acc

def reach_Mars(t, params):
    r = params[0:2]
    x = float(r[0])
    y = float(r[1])
    x_M = R_MARS*np.cos(ANGLE_V_MARS*t + OFFSET_MARS)
    y_M = R_MARS*np.sin(ANGLE_V_MARS*t + OFFSET_MARS)
    R2Mo = np.sqrt(x**2 + y**2) - R_MARS
    res = R2Mo 
    # print(res)
    # R2M = np.sqrt((x-x_M)**2 + (y-y_M)**2)
    # res = R2M 
    # geostationary_orbit = 3389500 + 17000000
    # R_MARS
    # print(R2M, geostationary_orbit)
    # if R2M <= geostationary_orbit:
        # res = 0
    return res

def model(t, params):
    r = params[0:2]
    v = params[2:]
    acceleration = (acc_g(t, r))
    # acceleration = (acc_g(t, r) + acc_dV(t, 0) + acc_dV(t, 1))
    # print(acc_dV(t, 0), acc_dV(t,1), t)
    return np.hstack((v, acceleration))


# t_end = 1.9*10**7
t_end = 3*10**7

reach_Mars.terminal = True
time_points = np.linspace(0, t_end, 10000)
sol = solve_ivp(model, [0, t_end], [R_EARTH + 35.786*10**6, 0, 0, V_EARTH + 3065 + dV1],
                t_eval=time_points,
                method='Radau',
                events=reach_Mars,
                max_step = 10000)
print(sol.t_events)

data = sol.y
x = data[0]
y = data[1]
t = sol.t

theta = np.arctan2(y,x)
R = np.sqrt(x**2 + y**2)
# plt.polar()


theta_E = ANGLE_V_EARTH*t
R_E = R_EARTH*np.ones(len(theta_E))
theta_M = ANGLE_V_MARS*t + OFFSET_MARS
R_M = R_MARS*np.ones(len(theta_M))


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
ax.plot(theta, R)
ax.plot(theta[-1], R[-1], 'o')
ax.plot(theta_E, R_E, color=[0,1,0])
ax.plot(theta_E[-1], R_E[-1], 'o', color=[0,1,0])
ax.plot(theta_M, R_M, color=[1,0,0])
ax.plot(theta_M[-1], R_M[-1], 'o', color=[1,0,0])
# plt.ylim(0, 1.5*R_EARTH)
# plt.plot(t,x)
# plt.plot(t,y)

x_E = R_EARTH*np.cos(theta_E)
y_E = R_EARTH*np.sin(theta_E)
x_M = R_MARS*np.cos(theta_M)
y_M = R_MARS*np.sin(theta_M)
# ax2.plot(t, x_E, 'g-')
# ax2.plot(t, x, 'k-')
ax2.plot(t, x_E, 'g-', t, y_E, 'g--')
ax2.plot(t, x_M, 'r-', t, y_M, 'r--')
ax2.plot(t, x, 'k-', t, y, 'k--')
# ax2.set_yscale('log')
ax2.grid()

ax3.plot(t, np.sqrt((x-x_M)**2 + (y-y_M)**2))
print('min', np.min(np.sqrt((x-x_M)**2 + (y-y_M)**2)))
ax3.set_yscale('log')
ax3.grid()

plt.show()