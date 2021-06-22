# Задание 11. Создайте модель процесса остывания стеклянного стакана с горячим кофе при 
# комнатных условиях. Постройте графики изменения температуры с учетом теплопроводности, 
# конвекции и испарения, а также оцените точность интегрирования в зависимости от схемы 
# интегрирования и величины шага интегрирования.

from matplotlib.colors import Colormap
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

AMBIENT_TEMP = 300 # Комнатная температура, К
HUMIDITY = 85 # Влажность, %
V_AIR = 0.1 # Скорость воздуха m/s

# Размеры сетки, в метрах
X_ROOM = 0.12 # m
Y_ROOM = 0.16 # m

# Температуропроводности веществ
T_COND_WATER = 0.143*10**(-6)
T_COND_GLASS = 3.4*10**(-7)
T_COND_AIR = 1.9*10**(-5)

L = 2260000 # Удельная теплота парообразования дж/кг
Cp = 4200 # Удельная теплоемкость, Дж/кг/К
rho = 1000 # Плотность воды кг/м3

# Давление насыщенного пара воды, мм.рт.ст
P_v_exp = np.array([4.585, 6.545, 9.212, 12.79, 17.54, 23.77, 31.84, 42.20, 55.37, 71.93, 92.59, 118.1])
P_v_exp_temps = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]) + 273

# Функция для расчета давления пара для заданой температуры
def get_Pv(T):
    return np.interp(T, P_v_exp_temps, P_v_exp)

coffee_start_temp = 350 # Начальная температура кофе, К

# Внутренние размеры жидкости
x_coffee = 0.04 #m
y_coffee = 0.1 #m

# Координаты стакана
x_glass1 = 0.04
x_glass2 = 0.05
y_glass1 = 0
y_glass2 = y_coffee

# Шаг координатной сетки
x_step = 0.005
y_step = 0.005


def get_laplacian(T):
    ''' Рассчет лапласиана для уравнения теплопроводности
    '''
    d2Tdx2 = np.zeros_like(T)
    d2Tdy2 = np.zeros_like(T)
    size_x = len(T[0])
    size_y = len(T[:,0])
    for j in range(0, size_y):
        for i in range(0, size_x):
            if i==0:
                d2Tdx2[j, i] = ( 2*T[j, i+1] - 2*T[j, i]) / x_step**2
            elif i==size_x-1:
                d2Tdx2[j, i] = ( -2*T[j, i] + 2*T[j, i-1] ) / x_step**2
            else:
                d2Tdx2[j, i] = ( T[j, i+1] - 2*T[j, i] + T[j, i-1] ) / x_step**2
            
            if j==0:
                d2Tdy2[j, i] = ( 2*T[j+1, i] - 2*T[j, i]) / y_step**2
            elif j==size_y-1:
                d2Tdy2[j, i] = ( -2*T[j, i] + 2*T[j-1, i] ) / y_step**2
            else:
                d2Tdy2[j, i] = ( T[j+1, i] - 2*T[j, i] + T[j-1, i] ) / y_step**2

    return d2Tdx2 + d2Tdy2


def get_evaporation_speed(T):
    # Скорость испарения, кг/с
    W = x_coffee*x_step * (0.022 + 0.0174 * V_AIR)*get_Pv(T)/1000*(1 - HUMIDITY/100) / 3600
    # W = np.pi*x_coffee**2 * (0.022 + 0.0174 * V_AIR)*get_Pv(T)*(1 - HUMIDITY/100) / 3600
    return W


def heat_equation(t, params):
    T = params[0:-1]
    T = T.reshape(initial_shape)
    m = params[-1]
    # Учет испарения
    dmdt = get_evaporation_speed(T[j_bnd, int(i_bnd/2)])
    dT = np.zeros_like(T)
    dT_evap = L*m/Cp/(rho*np.pi*x_coffee**2*y_step)
    dT[j_bnd,0:i_bnd] = dT_evap
    # Учет конвекции со стенкой
    lamb = 0.6
    for j in range(0, len(T[:,0])-1):
        if j < j_bnd:
            T[j, i_bnd] = (1 + x_step*120/lamb) * T[j, i_bnd-1] / 2

    dTdt = T_COND*get_laplacian(T) - dT

    return np.hstack((dTdt.flatten(), dmdt))


x = np.arange(0, X_ROOM, x_step)
y = np.arange(0, Y_ROOM, y_step)

i_bnd = None
j_bnd = None

X, Y = np.meshgrid(x, y)
start_temp_distr = np.zeros((len(y), len(x))) + AMBIENT_TEMP
T_COND = np.zeros_like(start_temp_distr)
# Задаем начальное распределение температуры, температуропроводности
for i, valx in enumerate(x):
    if not i_bnd:
        if valx>=x_coffee:
            i_bnd = i
            print(i_bnd)
    for j, valy in enumerate(y):
        if not j_bnd:
            if valy>=y_coffee:
                j_bnd = j
                print(j_bnd)
        if valx<x_coffee and valy<y_coffee:
            start_temp_distr[j, i] = coffee_start_temp
            T_COND[j, i] = T_COND_WATER
        elif valx>=x_glass1 and valx<=x_glass2 and valy>=y_glass1 and valy<=y_glass2:
            T_COND[j, i] = T_COND_GLASS
        else:
            T_COND[j, i] = T_COND_AIR

# Время окончания моделирования, с
t_end = 60
time_points = np.arange(0, t_end, 30)
initial_shape = start_temp_distr.shape
print("Init shape", initial_shape)
sol1 = solve_ivp(heat_equation, [0, t_end], np.hstack((start_temp_distr.flatten(), 0)),
                t_eval=time_points,
                method='Radau')
T = sol1.y
t = sol1.t
Ts = sol1.y[:,-1]
Ts = Ts[0:-1].reshape(initial_shape)


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X*100, Y*100, Ts-273, cmap='Reds',
                       linewidth=0, antialiased=True)

ax.set_ylabel('y, cm')
ax.set_xlabel('x, cm')
ax.set_zlabel('T, $^oC$')

plt.show()