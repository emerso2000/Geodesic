import sympy as sym
import numpy as np
from sympy import IndexedBase, Idx
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

        christoffel_symbols = []

        for i in range(4):
            for j in range(4):
                for k in range(4):
                    christoffel_symbols.append(0.5 * (sym.diff(self.g[i, j], coord_list[k]) +
                                                      sym.diff(self.g[i, k], coord_list[j]) -
                                                      sym.diff(self.g[j, k], coord_list[i])))

        christoffel_symbols = np.array(christoffel_symbols)
        christoffel_symbols = np.reshape(christoffel_symbols, (4, 4, 4))

        # print(christoffel_symbols)
        return christoffel_symbols

    def calc_geodesic_equations(self):
        """Calculate the geodesic equations for the Schwarzschild metric."""

        coord_list = [self.t, self.r, self.theta, self.phi]

        christoffel_symbols = self.calc_christoffel()

        geodesic_equations = []

        # Calculate the geodesic equations
        for i in range(4):
            geodesic_equations.append(0)
            for j in range(4):
                for k in range(4):
                    geodesic_equations[i] += christoffel_symbols[i, j, k] * sym.Derivative(
                        sym.Derivative(coord_list[k], coord_list[j]), coord_list[0])

        # print(geodesic_equations)
        return geodesic_equations

    # def calc_geodesic(self, initial_condition, step_size, num_steps):
    #     # coord_array = np.array([self.t, self.r, self.theta, self.phi])
    #     geodesic_equations = self.calc_geodesic_equations()
    #
    #     trajectory = sym.zeros(num_steps + 1, 4)
    #
    #     # Set the initial condition
    #     trajectory[0] = initial_condition
    #     for i in range(num_steps):
    #         # Check if the value of r has become equal to or less than the Schwarzschild radius
    #         # print(trajectory)
    #         if trajectory[0] <= self.r_s:
    #             print("Geodesic terminated at step", i, "due to r <= r_s")
    #             print(trajectory[i, 1])
    #             print(self.r_s)
    #             print(bool(trajectory[i, 1] <= self.r_s))
    #             break
    #         for j in range(4):
    #             trajectory[i + 1] = trajectory[i, j] + step_size * geodesic_equations[j].subs(
    #                 [(self.t, trajectory[i, 0]), (self.r, trajectory[i, 1]),
    #                  (self.theta, trajectory[i, 2]), (self.phi, trajectory[i, 3])])
    #
    #     #print(trajectory)
    #     return trajectory
    def calc_geodesic(self, initial_condition, step_size, num_steps):
        """Calculate the trajectory of a geodesic using Euler's method."""
        # Get the geodesic equations
        geodesic_equations = self.calc_geodesic_equations()

        # Initialize the trajectory array with the initial condition
        trajectory = sym.zeros(num_steps + 1, 4)
        trajectory[0] = initial_condition

        # Calculate the trajectory
        for i in range(num_steps):
            for j in range(4):
                trajectory[i + 1, j] = trajectory[i, j] + step_size * geodesic_equations[j].subs(
                    [(self.t, trajectory[i, 0]), (self.r, trajectory[i, 1]),
                     (self.theta, trajectory[i, 2]), (self.phi, trajectory[i, 3])])

        print(trajectory)
        return trajectory


initial_condition = sym.Matrix([0, 10, 0, 0])
schwarzschild = Schwarzschild(2, 6.67408e-11, 1.989e12)
schwarzschild.calc_christoffel()
schwarzschild.calc_geodesic_equations()
schwarzschild.calc_geodesic(initial_condition, 0.1, 1000)