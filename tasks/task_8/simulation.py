
import matplotlib.pyplot as plt
import numpy as np
from SolarSystem import System, Planets, Comet
from configs import configs

# Параметры кометы либо выбор заранее заданной конфигурации из configs.py,
# либо указать use_config = None для ручной настройки.
use_config = None
if use_config:
    cfg = configs[use_config]
    TheComet = Comet("Comet", cfg[0], cfg[1], cfg[2], cfg[3], cfg[4])
else:
    TheComet = Comet("Comet",
                     8*10**22,   # Масса, кг
                     1.5*10**12,  # Начальный радиус, м
                     0,  # Начальный угол, радианы
                     0,          # Начальная радиальная скорость, м/c
                     10000        # Начальная угловая скорость, м/с
                     )

### Параметры симуляции ###

# Время моделирования в секундах
t_end = 20*10**8
# Шаг времени в секундах  
dt = 2*10**4
# Временные отсчеты в которых посчитать значение координат
t = np.arange(0, t_end, dt)

# Моделируем движение кометы в Солнечной системе
# Метод интегрирования в функцие np.solve_ivp()
integration_method = 'Radau'
res = TheComet.evaluate_model([0, t_end], System, t_eval=t, method=integration_method)
res_time = res.t
comet_coords = res.y[0:2]
comet_vels = res.y[2:]
# Конвертируем координаты в полярные для дальнейшего построения графиков
comet_r, comet_theta = TheComet.ccdtp(comet_coords[0], comet_coords[1])
comet_vr, comet_vtheta = TheComet.csdtp(comet_coords[0], comet_coords[1], comet_vels[0], comet_vels[1])

### Далее код отвечает за отрисовку ###
fig, ax = plt.subplots(figsize=(8, 6), dpi=100, subplot_kw={'projection': 'polar'})
fig2, ax2 = plt.subplots(figsize=(8, 6), dpi=100)
fig3, ax3 = plt.subplots(figsize=(8, 6), dpi=100)

# Ресетим планеты в их изначальное положение
for Planet in Planets:
    Planet.reset()

# Расчет движения планет во времени (Нужно для для построения траекторий планет)
for i, val in enumerate(t):
    for Planet in Planets:
        r_next, theta_next = Planet.calculate_next_posistion(dt)

# Строим траектории планет и кометы в полярных координатах
for Planet in Planets:
    ax.plot(Planet.recorded_theta, Planet.recorded_r, '--', color=Planet.color, label=Planet.name)
    ax.plot(Planet.recorded_theta[-1], Planet.recorded_r[-1], 'o', color=Planet.color)

    ax2.plot(t, Planet.recorded_r[0:-1], label=Planet.name, color=Planet.color)

ax.plot(comet_theta, comet_r, color='k', label="Comet")
ax.plot(comet_theta[-1], comet_r[-1], 'o', color='k')
ax.plot(0, 0, 'o', color='#FFDF00', label="Sun")
ax.legend()

ax2.plot(res_time, comet_r, color='k', label="Comet")
ax2.grid(which='both')
ax2.legend(loc=3)
ax2.set_yscale('log')
ax2.set_xlabel("Время, с")
ax2.set_ylabel("Радиус траектории, м")

ax3.plot(res_time, comet_vr, label='$V_r$')
ax3.plot(res_time, comet_vtheta, label='$V_{\\theta}$')
ax3.set_xlabel("Время, с")
ax3.set_ylabel("Скорость, м/c")
ax3.grid(which='both')
ax3.legend(loc=3)

if use_config:
    fig.savefig('imgs_8/trj{}.png'.format(use_config), bbox_inches='tight')
    fig2.savefig('imgs_8/r{}.png'.format(use_config), bbox_inches='tight')
    fig3.savefig('imgs_8/v{}.png'.format(use_config), bbox_inches='tight')
else:
    fig.savefig('imgs_8/trj_custom.png', bbox_inches='tight')
    fig2.savefig('imgs_8/r_custom.png', bbox_inches='tight')
    fig3.savefig('imgs_8/v_custom.png', bbox_inches='tight')
plt.show()
