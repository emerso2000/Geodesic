import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

        # Calculate the inverse metric tensor
        g_inv = self.g.inv()

        christoffel_symbols = []

        for i in range(4):
            for j in range(4):
                for k in range(4):
                    sum = 0
                    for m in range(4):
                        sum += (0.5 * (sym.diff(self.g[i, j], coord_list[k]) +
                                       sym.diff(self.g[i, k], coord_list[j]) -
                                       sym.diff(self.g[j, k], coord_list[i]))
                                * g_inv[k, m])
                    christoffel_symbols.append(sum)

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

    # def calc_geodesic(self, initial_condition, step_size, num_steps):
    #     """Calculate the trajectory of a geodesic using Euler's method."""
    #     # Get the geodesic equations
    #     geodesic_equations = self.calc_geodesic_equations()
    #
    #     # Convert the symbolic equations into numerical functions
    #     num_geodesic_equations = [sym.lambdify((self.t, self.r, self.theta, self.phi), eq) for eq in geodesic_equations]
    #
    #     # Initialize the trajectory array with the initial condition
    #     trajectory = sym.zeros(num_steps + 1, 4)
    #     trajectory[0] = initial_condition
    #
    #     # Calculate the trajectory
    #     for i in range(num_steps):
    #         for j in range(4):
    #             # Substitute the values of the coordinates into the numerical geodesic equations
    #             # and evaluate the resulting equation
    #             trajectory[i + 1, j] = trajectory[i, j] + step_size * num_geodesic_equations[j](*trajectory[i])
    #
    #     return trajectory


initial_condition = sym.Matrix([0, 10, 0, 0])
schwarzschild = Schwarzschild(2, 6.67408e-11, 1.989e12)
print(schwarzschild.calc_christoffel())
print(schwarzschild.calc_geodesic_equations())
# schwarzschild.calc_geodesic(initial_condition, 0.1, 1000)
