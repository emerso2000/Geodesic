# Import libraries
import sympy as sym
import numpy as np
from sympy import IndexedBase, Idx
import matplotlib.pyplot as plt
import matplotlib.animation as animation


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

        # Save the metric tensor as an attribute of the Schwarzschild instance
        self.g = g

    def calc_christoffel(self):
        """Calculate the Christoffel symbols for the Schwarzschild metric."""

        # Create a list of the spacetime coordinates
        coord_list = [self.t, self.r, self.theta, self.phi]

        # Create a list of the Christoffel symbols
        christoffel_symbols = []

        # Calculate the Christoffel symbols
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    christoffel_symbols.append(0.5 * (sym.diff(self.g[i, j], coord_list[k]) +
                                                      sym.diff(self.g[i, k], coord_list[j]) -
                                                      sym.diff(self.g[j, k], coord_list[i])))

        # Convert the list of Christoffel symbols to a 4x4x4 array
        christoffel_symbols = np.array(christoffel_symbols)
        christoffel_symbols = np.reshape(christoffel_symbols, (4, 4, 4))

        print(christoffel_symbols)

    def calc_geodesic(self):
        pass
    def print_metric(self):
        """Print the metric tensor for the Schwarzschild metric."""
        print(self.g)


Schwarzschild = Schwarzschild(2, 6.67408e-11, 1.989e30)
Schwarzschild.print_metric()
Schwarzschild.calc_christoffel()
