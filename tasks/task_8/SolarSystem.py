# Задание 8. Постройте модель Солнечной системы. Рассчитайте параметры траектории кометы,
# попавшей в Солнечную систему извне. Постройте зависимости скорости и координаты кометы
# от времени при различных начальных параметрах, а также оцените точность интегрирования в
# зависимости от схемы интегрирования и величины шага интегрирования.
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

G = 6.67430*10**(-11)
M0 = 1.9885*10**30

class Sun():
    def __init__(self):
        self.mass = M0
        self.r = 0
        self.theta = 0
    
    def calculate_position_by_time(self, t):
        self.r = 0
        self.theta = 0
        return self.r, self.theta


class Comet():
    def __init__(self, name, mass, start_r, start_theta, v0_r, v0_theta):
        self.name = name
        self.mass = mass
        self.start_r = start_r
        self.start_theta = start_theta
        self.r = start_r
        self.theta = start_theta
        self.v_r = v0_r
        self.v_theta = v0_theta
        self.x, self.y = self.convert_coord_polar_to_decart(start_r, start_theta)
        self.v_x, self.v_y = self.convert_speed_polar_to_decart(v0_r, v0_theta, start_theta)
        self.recorded_theta = [start_theta]
        self.recorded_r = [start_r]
    
    def convert_coord_polar_to_decart(self, r, theta):
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return x, y

    def convert_speed_polar_to_decart(self, v_r, v_theta, theta):
        v_x = v_r*np.cos(theta) - v_theta*np.sin(theta)
        v_y = v_r*np.sin(theta) + v_theta*np.cos(theta)
        return v_x, v_y
    
    def convert_coord_decart_to_polar(self, x, y):
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return r, theta

    def calculate_distance_to_body(self, body):
        distance = np.sqrt( body.r**2 + self.r**2  - 2*body.r*self.r*np.cos(body.theta-self.theta))
        return distance
    
    def calculate_force_to_body(self, body):
        distance = self.calculate_distance_to_body(body)
        F = G*self.mass*body.mass/distance**2
        body_x, body_y = self.convert_coord_polar_to_decart(body.r, body.theta)
        F_x = F*(body_x - self.x)/distance
        F_y = F*(body_y - self.y)/distance
        return F_x, F_y

    def calculate_summary_acceleration(self, bodies):
        a_x_list, a_y_list = [], []
        for body in bodies:
            F_x, F_y = self.calculate_force_to_body(body)
            a_x_list.append(F_x/self.mass)
            a_y_list.append(F_y/self.mass)
        a_x = np.mean(a_x_list)
        a_y = np.mean(a_y_list)
        return a_x, a_y

    def model_func(self, t, data_vec):
        #  data_vec = [x,y, vx,vy]
        coords = data_vec[0:2]
        vel = data_vec[2:]

        self.x, self.y = coords[0], coords[1]
        self.r = np.sqrt(coords[0]**2 + coords[1]**2)
        self.theta = np.arctan2(coords[1], coords[0])

        for body in self.bodies:
            body.calculate_position_by_time(t)
        ax, ay = self.calculate_summary_acceleration(self.bodies)
        dvdt = [ax, ay]
        # returns r', v'
        return np.hstack((vel, dvdt))

    def evaluate_model(self, t_end, bodies, t_eval=None, method='RK45'):
        self.bodies = bodies
        res = solve_ivp(self.model_func, t_end,
                        y0=[self.x, self.y, self.v_x, self.v_y],
                        t_eval=t_eval, method=method)
        return res


class CelestialBody():
    def __init__(self, name, mass, a, eccentricity,
                 offset, start_theta, plot_color):
        ''' Основной класс для описания орбит планет вокруг Солнца.
            Включает методы для рассчета движения во времени по заданной траектории.
        '''
        self.name = name
        self.mass = mass 
        self.a = a  # Большая полуось
        self.ecc = eccentricity  # Эксцентриситет
        self.angle_offset = np.pi*offset/180 # Долгота восходящего узла
        self.start_theta = start_theta # Начальное положение
        self.p = a*(1-self.ecc**2)
        self.mu = G*M0
        self.color = plot_color

        self.const1 = self.a*(1-self.ecc**2)
        self.const2 = np.sqrt(self.mu*self.a*(1-self.ecc**2))
        self.get_r(self.start_theta)
        self.theta = start_theta

        self.recorded_theta = [self.start_theta]
        self.recorded_r = [self.r]
    
    def calculate_position_by_time(self, t, points=100):
        dt = t/points
        for i in np.linspace(0, t, points):
            self.calculate_next_posistion(dt)
        r, theta = self.recorded_r[-1], self.recorded_theta[-1]
        self.reset
        self.r, self.theta = r, theta
        return r, theta

    def get_r(self, theta):
        ''' Get radius value for current angle theta.
        '''
        self.r = self.const1/(1+self.ecc*np.cos(theta + self.angle_offset))
        return self.r

    def get_d_theta(self, dt):
        '''  Get angle increment d_theta for current radius value
             and time increment dt.
        '''
        d_theta = 1/self.r**2*self.const2*dt
        return d_theta

    def calculate_next_posistion(self, dt):
        d_theta = self.get_d_theta(dt)
        next_theta = self.recorded_theta[-1] + d_theta
        next_r = self.get_r(next_theta)
        self.recorded_theta.append(next_theta)
        self.recorded_r.append(self.r)
        return next_r, next_theta

    def reset(self):
        self.get_r(self.start_theta)
        self.recorded_theta = [self.start_theta]
        self.recorded_r = [self.r]

TheSun = Sun()

Mercury = CelestialBody("Mercury", 3.33022*10**23, 57909227000, 0.20563593, 48.33167, -0.2*np.pi, [1,0.5,0])
Venus = CelestialBody("Venus", 4.8675*10**24, 108208930000, 0.0068, 76.67069, 0.1*np.pi, [1,0,0.4])
Earth = CelestialBody("Earth", 5.9726*10**24, 149598261000, 0.01671123, 348.73936, 0.8*np.pi, [0,1,0.3])
Mars = CelestialBody("Mars", 6.4171*10**23, 2.2794382*10**8*1000, 0.0933941, 49.57854, 1.15*np.pi, [1,0,0])
Jupiter = CelestialBody("Jupiter", 1.8986*10**27, 7.785472*10**8*1000, 0.048775, 100.55615, -0.02*np.pi, [0.9,0.7,0.3])
Saturn = CelestialBody("Saturn", 5.6846*10**26, 1429394069000, 0.055723219, 113.642, -0.08*np.pi, [0.6, 0.1, 0.3])
Uranus = CelestialBody("Uranus", 8.6813*10**25, 2876679082000, 0.044, 73.9898, 0.06*np.pi, [0,0.5,0.9])
Neptune = CelestialBody("Neptune", 1.0243*10**26, 4503443661000, 0.011214269, 131.794, 0.3*np.pi, [0,0.1,1])
Planets = [Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
System = [TheSun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
