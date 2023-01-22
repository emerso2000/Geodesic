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
        self.x = 0
        self.y = 0.02
        self.z = 0.01
        
        self.angle = angle_initial
        self.counter = counter_initial
        self.points = []
        
        self.a = 10.0
        self.b = 28.0
        self.c = 8.0 / 3.0

    def lorenz_equations(self):
        dx = 0
        dy = 0
        dz = 0

        # glVertex3d(self.x, self.y, self.z)

        self.x = self.x + dx
        self.y = self.y + dy
        self.z = self.z + dz

        self.counter += 1

        self.points.append((self.x, self.y, self.z))
        

integrator = LorenzAttractor()
print(integrator.points)
