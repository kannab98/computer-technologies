
import matplotlib.pyplot as plt
import numpy as np
from SolarSystemOriginal import System, TheSun, Planets, TheComet

# dt = 1
dt = 100000
# t = np.arange(0, 10, dt)
t = np.arange(0, 1050050000, dt)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# Строим траектории планет
# theta = np.linspace(0, 2*np.pi, 1000)
# for Planet in Planets:
#     # ax.plot(Planet.start_theta, Planet.r, 'o', color=Planet.color)
#     R = Planet.get_r(theta)
#     ax.plot(theta, R, color=Planet.color, alpha=0.3)
#     # ax.plot(theta, R, color=Planet.color, label=Planet.name, alpha=0.3)
#     Planet.reset()

# Симуляция движения планет во времени
for i, val in enumerate(t):

    TheComet.calculate_next_position(dt)
    # a_x, a_y = TheComet.calculate_summary_acceleration([TheSun])
    a_x, a_y = TheComet.calculate_summary_acceleration(System)

    new_v_x = TheComet.v_x + a_x*dt
    new_v_y = TheComet.v_y + a_y*dt
    TheComet.v_x = new_v_x
    TheComet.v_y = new_v_y

    for Planet in Planets:

        r_next, theta_next = Planet.calculate_next_posistion(dt)

for Planet in Planets:
    ax.plot(Planet.recorded_theta, Planet.recorded_r, '--', color=Planet.color, label=Planet.name)
    ax.plot(Planet.recorded_theta[-1], Planet.recorded_r[-1], 'o', color=Planet.color)

ax.plot(TheComet.recorded_theta, TheComet.recorded_r, color='k', label="Comet")
ax.plot(TheComet.recorded_theta[-1], TheComet.recorded_r[-1], 'o', color='k')
plt.legend()
plt.show()
