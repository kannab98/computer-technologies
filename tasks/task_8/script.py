
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from SolarSystem import System, TheSun, Planets, TheComet

# dt = 1
dt = 1000000
t_end = 15050050000
# 3850050000
# t = np.arange(0, 10, dt)
t = np.arange(0, t_end, dt)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()


res = TheComet.evaluate_model([0, t_end], System, method='RK23')
# res = TheComet.evaluate_model([0, t_end], [TheSun], t_eval=t, method='RK23')
time = res.t
coords = res.y[0:2]
comet_r, comet_theta = TheComet.convert_coord_decart_to_polar(coords[0], coords[1]) 
vels = res.y[2:]

ax3.plot(time, coords[0], label='x')
ax3.plot(time, coords[1], label='y')

ax2.plot(time, vels[0], label='Vx')
ax2.plot(time, vels[1], label='Vy')

for Planet in Planets:
    Planet.reset()
# Симуляция движения планет во времени
for i, val in enumerate(t):
    for Planet in Planets:
        r_next, theta_next = Planet.calculate_next_posistion(dt)

for Planet in Planets:
    ax.plot(Planet.recorded_theta, Planet.recorded_r, '--', color=Planet.color, label=Planet.name)
    ax.plot(Planet.recorded_theta[-1], Planet.recorded_r[-1], 'o', color=Planet.color)
    # ax.plot(Planet.last_theta, Planet.last_r, '*', color='k')

ax.plot(comet_theta, comet_r, color='k', label="Comet")
ax.plot(comet_theta[-1], comet_r[-1], 'o', color='k')
ax.plot(0, 0, 'o', color='#FFDF00', label="Sun")
plt.legend()
plt.grid()
plt.show()