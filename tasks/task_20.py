import numpy as np
from numpy import random


class Tor(object):
    def __init__(self, a: tuple) -> None:
        # R is the distance from the center of the tube to the center of the torus,
        self.R = a[0]
        # r is the radius of the tube.
        self.r = a[1]

    def mask(self, r: float) -> bool:

        r0 = self.r
        R0 = self.R

        # XY-axis radius projection
        rho2 = np.sum(r[0:2]**2, axis=0)
        # Radius absolute value
        R2 = rho2 + r[2]**2
        # Surface equation
        # https://en.wikipedia.org/wiki/Torus
        mask = (R2 + R0**2 - r0**2)**2 - 4*R0**2*(rho2) <= 0
        return mask

    @property
    def volume(self) -> float:
        # Analytical expression for volume
        return 2*np.pi**2*self.R*self.r**2

    def grid(self, num: tuple, dtype="linear"):
        R = self.R
        r = self.r

        if dtype == "linear":
            return cartesian_meshgrid([-R-r, -R-r, -r], [R+r, R+r, +r], num)
        elif dtype == "random":
            return cartesian_random([-R-r, -R-r, -r], [R+r, R+r, +r], np.prod(num))


def monte_carlo(mask: bool, r: float):
    axis = 1
    # Volume of an N-dimensional cube describing the area of integration
    area = np.prod(r.max(axis) - r.min(axis), dtype=np.float64)
    return np.mean(mask) * area


def rect(mask: bool, dr: float):
    # Volume of elementary N-dimensional cubes dV inside area of integration
    dV = np.prod(dr, axis=0, where=mask)
    return np.sum(dV, where=mask)


def cartesian_meshgrid(start, stop, num):
    x = np.linspace(start[0], stop[0], num[0])
    y = np.linspace(start[1], stop[1], num[1])
    z = np.linspace(start[2], stop[2], num[2])

    dx = np.zeros_like(x)
    dy = np.zeros_like(y)
    dz = np.zeros_like(z)

    dx[1:] = np.diff(x)
    dy[1:] = np.diff(y)
    dz[1:] = np.diff(z)

    x, y, z = np.meshgrid(x, y, z)
    dx, dy, dz = np.meshgrid(dx, dy, dz)

    return np.array([x, y, z]), np.array([dx, dy, dz])


def cartesian_random(start, stop, size):
    x = random.uniform(start[0], stop[0], size)
    y = random.uniform(start[1], stop[1], size)
    z = random.uniform(start[2], stop[2], size)
    return np.array([x, y, z])
