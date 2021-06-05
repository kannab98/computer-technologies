
import matplotlib.pyplot as plt
import numpy as np
from SolarSystem import System, Planets, Comet
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

TheComet = Comet("Comet",
                 1*10**11,   # Масса
                 50*10**11,  # Начальный радиус
                 3*np.pi/4,  # Начальный угол
                 0,          # Начальная радиальная скорость
                 1500        # Начальная угловая скорость
                 )

# Параметры симуляции
integration_method = 'RK23'
t_end = 1*10**9  
dt = 2*10**6
t = np.arange(0, t_end, dt)


# Моделируем модель кометы в Солнечной системе
res = TheComet.evaluate_model(
    [0, t_end], System, t_eval=t, method=integration_method)
time = res.t
coords = res.y[0:2]
vels = res.y[2:]
# Конвертируем координаты в полярные для дальнейшего построения графиков
comet_r, comet_theta = TheComet.convert_coord_decart_to_polar(
    coords[0], coords[1])

# Ресетим планеты в их изначальное положение
for Planet in Planets:
    Planet.reset()

# Расчет движения планет во времени (Нужно для для построения траекторий планет)
for i, val in enumerate(t):
    for Planet in Planets:
        r_next, theta_next = Planet.calculate_next_posistion(dt)

# Строим траектории планет и кометы в полярных координатах
for Planet in Planets:
    ax.plot(Planet.recorded_theta, Planet.recorded_r,
            '--', color=Planet.color, label=Planet.name)
    ax.plot(Planet.recorded_theta[-1],
            Planet.recorded_r[-1], 'o', color=Planet.color)
# ax.plot(comet_theta, comet_r, color='k', label="Comet")
# ax.plot(comet_theta[-1], comet_r[-1], 'o', color='k')
ax.plot(0, 0, 'o', color='#FFDF00', label="Sun")

plt.legend()
plt.show()