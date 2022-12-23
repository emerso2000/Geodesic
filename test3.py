import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import math


class Schwarzschild:
    def __init__(self, r_s, G, M):
        self.r_s = r_s  # Schwarzschild radius
        self.G = G  # Gravitational constant
        self.M = M  # Mass of central object

        # Define the spacetime coordinates (t, r, theta, phi)
        self.t, self.r, self.theta, self.phi = sym.symbols('t r theta phi')

        # Define the metric tensor components in terms of the Schwarzschild radius, gravitational constant,
        # and mass of the central object
        self.g_tt = -(1 - self.r_s / self.r)
        self.g_rr = 1 / (1 - self.r_s / self.r)
        self.g_theta_theta = self.r ** 2
        self.g_phi_phi = (self.r * sym.sin(self.theta)) ** 2

        # Define the metric tensor as a 4x4 matrix using the components defined above
        g = sym.Matrix([[self.g_tt, 0, 0, 0],
                        [0, self.g_rr, 0, 0],
                        [0, 0, self.g_theta_theta, 0],
                        [0, 0, 0, self.g_phi_phi]])

        self.g = g

    def calc_christoffel(self):
        """Calculate the Christoffel symbols for the Schwarzschild metric."""

        coord_list = [self.t, self.r, self.theta, self.phi]
        g_inv = self.g.inv()
        christoffel_symbols = []

        for i in range(4):
            for j in range(4):
                for k in range(4):
                    total = 0
                    for m in range(4):
                        total += (0.5 * (sym.diff(self.g[i, m], coord_list[k]) +
                                         sym.diff(self.g[i, m], coord_list[j]) -
                                         sym.diff(self.g[i, m], coord_list[i]))
                                  * g_inv[m, k])
                    christoffel_symbols.append(total)

        christoffel_symbols = np.array(christoffel_symbols)
        christoffel_symbols = np.reshape(christoffel_symbols, (4, 4, 4))

        return christoffel_symbols

    def calc_geodesic_equations(self):
        """Calculate the geodesic equations for the Schwarzschild metric."""
        christoffel_symbols = self.calc_christoffel()

        u_t, u_r, u_theta, u_phi = sym.symbols('u_t u_r u_theta u_phi')

        u = sym.Matrix([u_t, u_r, u_theta, u_phi])

        a = sym.zeros(4, 1)

        for mu in range(4):
            for alpha in range(4):
                for beta in range(4):
                    a[mu] += (-christoffel_symbols[mu][alpha][beta] * u[alpha] * u[beta] +
                              christoffel_symbols[0][alpha][beta] * u[alpha] * u[beta] * u[mu])

        return a.tolist()


class Point:
    def __init__(self, t, r, theta, phi):
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi


class Body:
    def __init__(self, location, mass, velocity, name=""):
        self.location = location
        self.mass = mass
        self.velocity = velocity
        self.name = name


class EulerIntegrator:
    def __init__(self, timestep, bodies, r_s, G, M):
        self.time_step = timestep
        self.bodies = bodies
        self.r_s = r_s
        self.G = G
        self.M = M
        self.geodesic_equations = Schwarzschild(self.r_s, self.G, self.M).calc_geodesic_equations()

    def calc_acceleration(self, body):
        """Calculate the acceleration of a body using Euler's method."""

        pass


integrator = EulerIntegrator(0.01, [], 1, 1, 1)
print(integrator.calc_acceleration(Body(Point(1, 1, 1, 1), 1, 1)))
