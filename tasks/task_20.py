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

def cube(a=2):
    x = np.linspace(0, a)
    y = np.linspace(0, a)
    z = np.linspace(0, a)
    x,y,z = np.meshgrid(x,y,z)
    return x,y,z

def monte_carlo(func: np.ndarray, x: np.ndarray, y: np.ndarray):
    V = np.sum(func)
    return V

def crect(func: np.ndarray, x: np.ndarray, y: np.ndarray):
    V = 0
    for i in range(1, func.shape[0]):
        for j in range(1, func.shape[1]):
            for k in range(1, func.shape[2]):
                V += (func[i,j,k] - func[i-1, j-1, k-1]) * \
                        (x[i,j,k] - x[i-1, j-1, k-1]) * \
                        (y[i,j,k] - y[i-1, j-1, k-1])
    return V 

a = 3
s = 50
x = np.linspace(0, a, s)
y = np.linspace(0, a, s)
z = np.linspace(0, a, s)
x,y,z = np.meshgrid(x,y,z)

V = np.sum(z)
V0 = crect(z,x,y)
print(V, V0)

