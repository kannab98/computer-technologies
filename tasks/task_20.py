import numpy as np

def tor(r=1, R=2):
    """
    Генерация тороидальной поверхности
    """
    el = np.linspace(-np.pi, np.pi, 32, endpoint=True)
    az = np.linspace(0, 2*np.pi, 32, endpoint=True)

    el, az = np.meshgrid(az, el)


    x = (R + r*np.cos(el)) * np.cos(az) 
    y = (R + r*np.cos(el)) * np.sin(az) 
    z = r * np.sin(el)

    return x, y, z

def monte_carlo():
    pass

def crect():
    pass