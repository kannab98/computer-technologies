# Задание 9. Рассчитайте эффективную траекторию ракеты, предназначенной для наименее
# затратного по топливу запуска с Земли искусственного спутника Марса. Постройте зависимости
# скорости и координаты ракеты от времени, а также оцените точность интегрирования в
# зависимости от схемы интегрирования и величины шага интегрирования

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Константы 
G = 6.67430*10**(-11)
M_SUN = 1.989*10**30  # kg

# Параметры Марса и его орбиты
M_MARS = 6.39*10**23  # кг
R_MARS = 2.26655*10**11  # Радиус орбиты, м
V_MARS = 24130  # Орбитальная скорость, м/с
ANGLE_V_MARS = V_MARS/R_MARS  # Угловая скорость, рад/с
OFFSET_MARS = 0.29*np.pi
GSO_MARS = 3389500 + 17000000 # Радиус ГСО от центра, м

# Параметры Земли и ее орбиты
M_EARTH = 5.972*10**24  # кг
R_EARTH = 150*10**9  # Радиус орбиты, м
V_EARTH = 29783  # Орбитальная скорость, м/с
ANGLE_V_EARTH = V_EARTH/R_EARTH  # Угловая скорость, рад/с
GSO_EARTH = 6371000 + 35.786*10**6 # Радиус ГСО от центра, м

M_ROCKET = 3*10**6 # Масса ракеты

# Теоретические величины импульсов 
Vp = np.sqrt((V_MARS**2 + V_EARTH**2)/2)
dV1 = V_EARTH*(V_EARTH/Vp - 1)
dV2 = V_MARS*(1 - V_MARS/Vp)

def acc_g(t, r, offset_mars):
    ''' Рассчет суммарного ускорения от трех тел - 
        Земли, Марса и Солнца.
    '''
    x, y = r
    x_E = R_EARTH*np.cos(ANGLE_V_EARTH*t)
    y_E = R_EARTH*np.sin(ANGLE_V_EARTH*t)
    x_M = R_MARS*np.cos(ANGLE_V_MARS*t + offset_mars)
    y_M = R_MARS*np.sin(ANGLE_V_MARS*t + offset_mars)

    # Расстояния до тел 
    R2E = np.sqrt((x-x_E)**2 + (y-y_E)**2)
    R2M = np.sqrt((x-x_M)**2 + (y-y_M)**2)
    R2S = np.sqrt(x**2 + y**2)
    
    # Ускорения от каждого тела 
    acc_E = M_EARTH/R2E**3 * np.array([x_E - x, y_E - y])
    acc_M = M_MARS/R2M**3 * np.array([x_M - x, y_M - y])
    acc_S = -M_SUN/R2S**3 * np.array([x, y])

    acc = G*(acc_E + acc_M + acc_S)
    return acc


def model(t, params, offset_mars):
    ''' Модель для моделирования - рассчет
        дифференциального уравнения движения
    '''
    r = params[0:2]
    v = params[2:]
    acceleration = (acc_g(t, r, offset_mars))
    return np.hstack((v, acceleration))


def find_optimal_parameters(params):
    ''' Функция для минимизации возвращаемого значения.
        Оптимальными будут те параметры, при которых 
        расстояние до Марса будет меньше величины ГСО, радиальная скорость
        в максимуме Орбиты будет минимальна, а также максимальная величины
        радиуса орбиты ракеты будет максимально близка к орбите Марса.
    '''
    dv, offset_mars = params
    # Расчет траекотрии при заданных параметрах
    sol = solve_ivp(model, [0, t_end],
                    [R_EARTH + GSO_EARTH, 0, 0, V_EARTH + 3065 + dv],
                    args=(offset_mars,),
                    t_eval=time_points,
                    method='Radau'
                    )
    data, t = sol.y, sol.t
    x, y, vx, vy = data
    theta = np.arctan2(y, x)
    R = np.sqrt(x**2 + y**2)
    vr = vx*np.cos(theta) + vy*np.sin(theta)
    
    theta_M = ANGLE_V_MARS*t + offset_mars
    x_M = R_MARS*np.cos(theta_M)
    y_M = R_MARS*np.sin(theta_M)
    # Расстояние до Марса в каждый момент времени
    d2M = np.sqrt((x-x_M)**2 + (y-y_M)**2)

    min_d2M = np.min(d2M)
    ind_min = np.argmin(d2M)
    vr_at_min = vr[ind_min] # Радиальная скорость при максимальном сближении с Марсом
    
    max_R = np.max(R)
    d2MO = np.abs(max_R - R_MARS) # Расстояние от максимальной точки орбиты ракеты до орбиты Марса
    min_d2M -= GSO_MARS # Минимальное расстояние до ГСО Марса

    # Вывод значений при оптимизации
    # print('Imp:', dv, '\tM off:', offset_mars, '\tmin vr:',
        #   vr_at_min, '\tmin d2M:', min_d2M, ' m')
    return np.abs(min_d2M) + vr_at_min**2 + d2MO

# Первичное время моделирования
t_end = 3*10**7
dV_initial = 2400 # Изначальное предположение первого импульса, м/с
time_points = np.linspace(0, t_end, 10000)

# Моделирование с поиском оптимальной величины начального импульса и положения Марса
print("Ищем оптимальные параметры запуска...")
res = minimize(find_optimal_parameters, [dV_initial, OFFSET_MARS],
               options={'maxiter':200}, method='Nelder-Mead')

print(res.message)
opt_start_impulse, opt_offset_mars = res.x # Найденные оптимальные значения импульса и положения Марса
print("Теоретическое значение импульса:", dV1, 'м/с')
print("Изначальное значение импульса:", dV_initial, 'м/с')
print("Оптимальное значение импульса:", opt_start_impulse, 'м/с')
print("Изначальное положение Марса:", OFFSET_MARS, 'радиан')
print("Оптимальное положение Марса:", opt_offset_mars, 'радиан')

# Рассчет траектории с оптимальным значением начального импульса скорости
start_values = [R_EARTH + GSO_EARTH, 0, 0, V_EARTH + 3065 + opt_start_impulse]
sol1 = solve_ivp(model, [0, t_end], start_values, args=(opt_offset_mars,),
                t_eval=time_points, method='Radau')
data1, t1 = sol1.y, sol1.t
x1, y1, vx1, vy1 = data1
theta = np.arctan2(y1, x1)

theta_M = ANGLE_V_MARS*t1 + opt_offset_mars
x_M = R_MARS*np.cos(theta_M)
y_M = R_MARS*np.sin(theta_M)

D2M = np.sqrt((x1-x_M)**2 + (y1-y_M)**2)
min_D2M = np.min(D2M) # Минимальное расстояние до Марса
min_D2M_ind = np.argmin(D2M)
print('Расстояние до Марса при максимальном сближении:', min_D2M/1000,
      'км, Расстояние ГСО Марса', GSO_MARS/1000, 'км')
ETA_M = t1[min_D2M_ind]  # Время максимального сближения

# На данном этапе нам интересны только время и координаты до сближения
list_of_data = [t1, x1, y1, vx1, vy1, theta, theta_M, x_M, y_M]
for i in range(0, len(list_of_data)):
    list_of_data[i] = list_of_data[i][0:min_D2M_ind+1]
t1, x1, y1, vx1, vy1, theta, theta_M, x_M, y_M = list_of_data

vr1 = vx1*np.cos(theta) + vy1*np.sin(theta)
vr_at_min = vr1[-1] # Радиальная скорость при максимальном сближении
print("Время до максимального сближения:", ETA_M/86400, 'дней')
print('Радиальная скорость при максимальном сближении:', vr_at_min, 'м/c')

# Расчет продолжения траектории 
print("Расчет продолжения тракетории")
print("Величина второго импульса:", dV2, "м/c")
extra_time = 0.5*10**7 # Дополнительное время для рассчета, с
time_points2 = np.linspace(ETA_M, ETA_M + extra_time, 5000)
start_values2 = [x1[-1], y1[-1], vx1[-1], vy1[-1] - dV2]
sol2 = solve_ivp(model, [ETA_M, ETA_M + extra_time], start_values2, args=(opt_offset_mars,),
                t_eval=time_points2, method='Radau')
data2, t2 = sol2.y, sol2.t
x2, y2, vx2, vy2 = data2

# Временный отсчеты
t = np.concatenate((t1, t2))

# Координаты ракеты
x = np.concatenate((x1, x2))
y = np.concatenate((y1, y2))
vx = np.concatenate((vx1, vx2))
vy = np.concatenate((vy1, vy2))
theta = np.arctan2(y, x)
R = np.sqrt(x**2 + y**2)

# Координаты Земли
theta_E = ANGLE_V_EARTH*t
R_E = R_EARTH*np.ones(len(theta_E))
x_E = R_EARTH*np.cos(theta_E)
y_E = R_EARTH*np.sin(theta_E)

# Координаты Марса
theta_M = ANGLE_V_MARS*t + opt_offset_mars
R_M = R_MARS*np.ones(len(theta_M))
x_M = R_MARS*np.cos(theta_M)
y_M = R_MARS*np.sin(theta_M)

#########    Блок отрисовки    ####### 
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

# Траектории движения
ax.plot(theta, R, label='Ракета')
ax.plot(theta[-1], R[-1], 'o')
ax.plot(theta_E, R_E, color=[0, 1, 0], label='Земля')
ax.plot(theta_E[-1], R_E[-1], 'o', color=[0, 1, 0])
ax.plot(theta_M, R_M, color=[1, 0, 0], label='Марс')
ax.plot(theta_M[-1], R_M[-1], 'o', color=[1, 0, 0])
ax.legend()

# Координаты
ax2.plot(t, x_E, 'g-', label='X Земля')
ax2.plot(t, y_E, 'g--', label='Y Земля')
ax2.plot(t, x_M, 'r-', label='X Марс')
ax2.plot(t, y_M, 'r--', label='Y Марс')
ax2.plot(t, x, 'k-', label='X Ракета')
ax2.plot(t, y, 'k--', label='Y Ракета')
ax2.axvline(ETA_M, label='Прибытие на Марс')
ax2.set_xlabel('Время, с')
ax2.set_ylabel('Координата, м')
ax2.legend()
ax2.grid()

# Скорость
ax3.plot(t, vx, 'k-', label='V_x')
ax3.plot(t, vy, 'k--', label='V_y')
ax3.axvline(ETA_M, label='Прибытие на Марс')
ax3.set_xlabel('Время, с')
ax3.set_ylabel('Скорость, м/с')
ax3.legend()
ax3.grid()

plt.show()
