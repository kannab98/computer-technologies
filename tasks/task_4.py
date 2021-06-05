import numpy as np
import matplotlib.pyplot as plt

# Решение систем дифференциальных уравнений
from scipy.integrate import solve_ivp

# Задание 4. Разработайте алгоритм решения задачи для модели, описывающей вынужденные
# колебания системы, включающей три тела, соединенные пружинами, при наличии силы вязкого
# сопротивления. Постройте осциллограммы и фазовые портреты колебаний, а также оцените
# точность интегрирования в зависимости от схемы интегрирования и величины шага
# интегрирования.


sol = solve_ivp(predator_prey, t_span=[0, 500], y0=[10, 2], args=(p,"McA"))

