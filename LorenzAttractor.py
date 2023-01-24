import pygame
from pygame.locals import *
import time
import math
from math import *
from OpenGL.GL import *
from OpenGL.GLU import *
import random

x_initial = 0
y_initial = 0.02
z_initial = 0.01

angle_initial = 0.01
display = (720, 720)

counter_initial = 0


class LorenzAttractor:
    def __init__(self):
        self.x = 0.01
        self.y = 0
        self.z = 0

        self.angle = 0.01
        self.counter = 0
        self.points = []

        self.a = 10.0
        self.b = 28.0
        self.c = 8.0 / 3.0

        self.display = (720, 720)

    def lorenz_equations(self):
        dx = 0
        dy = 0
        dz = 0

        dt = 0.01

        while self.counter < 250:
            glColor3d(1, 1, 1)
            dx = self.a * (self.y - self.x) * dt
            dy = (self.x * (self.b - self.z) - self.y) * dt
            dz = (self.x * self.y - self.c * self.z) * dt

            glVertex3d(self.x, self.y, self.z)

            self.x += dx
            self.y += dy
            self.z += dz

            self.counter += 1

            self.points.append((self.x, self.y, self.z))

        return self.points


integrator = LorenzAttractor()
points = integrator.lorenz_equations()

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
plt.title("Lorenz Attractor")
for x, y, z in points:
    ax.scatter3D(x, y, z)
plt.show()
