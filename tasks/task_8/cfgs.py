# Перехват Венерой
TheComet = Comet("Comet",
                 7.3477*10**22,   # Масса, кг
                 152458261000,  # Начальный радиус, м
                 0.8*np.pi,  # Начальный угол, радианы
                 0,          # Начальная радиальная скорость, м/c
                 29800        # Начальная угловая скорость, рад/с
                 )

### Параметры симуляции ###

# Метод интегрирования в функцие np.solve_ivp()
integration_method = 'Radau'
# Время моделирования в секундах
t_end = 1.354*10**7
# t_end = 3.154*10**7  # Земной год
# Шаг времени в секундах  
dt = 2*10**4
# Временные отсчеты в которых посчитать значение координат
t = np.arange(0, t_end, dt)


# Типа луна
TheComet = Comet("Comet",
                 7.3477*10**22,   # Масса, кг
                 150537211000,  # Начальный радиус, м
                 0.8*np.pi,  # Начальный угол, радианы
                 0,          # Начальная радиальная скорость, м/c
                 30000        # Начальная угловая скорость, рад/с
                 )

### Параметры симуляции ###

# Метод интегрирования в функцие np.solve_ivp()
integration_method = 'Radau'
# Время моделирования в секундах
t_end = 0.98*10**7
# t_end = 3.154*10**7  # Земной год
# Шаг времени в секундах  
dt = 2*10**4


# yonk
TheComet = Comet("Comet",
                 7.3477*10**22,   # Масса, кг
                 1550537211000,  # Начальный радиус, м
                 0.8*np.pi,  # Начальный угол, радианы
                 0,          # Начальная радиальная скорость, м/c
                 10000        # Начальная угловая скорость, м/с
                 )

### Параметры симуляции ###

# Метод интегрирования в функцие np.solve_ivp()
integration_method = 'Radau'
# Время моделирования в секундах
t_end = 20.98*10**8
