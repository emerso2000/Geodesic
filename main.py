import sympy as sp, numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Schwarzchild:
    def __init__(self):
        self.t, self.r, self.θ, self.ϕ, self.M = sp.symbols("t r θ ϕ M")
        self.vars = [self.t, self.r, self.θ, self.ϕ]

    def metric_tensor(self):
        g_co = sp.Matrix(
            [
                [-(1 - 2 * self.M / self.r), 0, 0, 0],
                [0, 1 / (1 - 2 * self.M / self.r), 0, 0],
                [0, 0, self.r ** 2, 0],
                [0, 0, 0, self.r ** 2 * sp.sin(self.θ) ** 2],
            ]
        )

        g_contra = g_co.inv()
        dg_co_by = sp.Array([[[sp.diff(g_co[i, j], self.vars[k]) for k in range(4)] for j in range(4)] for i in range(4)])

        return g_co, g_contra, dg_co_by

    def basic_sum(self, it):
        it = iter(it)
        cur = next(it)
        for el in it:
            cur += el
        return cur

    def christoffel_symbols(self):
        g_co, g_contra, dg_co_by = self.metric_tensor()
        Γ = sp.Array(
            [
                [
                    [
                        self.basic_sum(
                            g_contra[i, m] * (dg_co_by[m, k, l] + dg_co_by[m, l, k] - dg_co_by[k, l, m]) for m in
                            range(4)
                        )
                        / 2
                        for l in range(4)
                    ]
                    for k in range(4)
                ]
                for i in range(4)
            ]
        )
        return Γ

    def geodesic_equation(self):
        s = sp.symbols("s")
        Γ = self.christoffel_symbols()

        (*dvars,) = dt, dr, dθ, dϕ = sp.symbols("dt dr dθ dϕ")  # "speeds" (derivatives by ds) of each coord
        accel = sp.Array(
            [
                -self.basic_sum(Γ[mu, alpha, beta] * (dvars[alpha] * dvars[beta]) for alpha in range(4) for beta in range(4))
                for mu in range(4)
            ],
            shape=(4, 1),
        )
        return accel

    def calc_time_component(self, positions, velocity):
        r = positions[1]
        sint = positions[2]
        r_s = 1

        rVel = velocity[0]
        thetaVel = velocity[1]
        phiVel = velocity[2]

        epsilon = 0

        timeVel = math.sqrt((epsilon + (1 / (1 - (r_s / r))) * rVel ** 2) + (
                (r ** 2) * ((thetaVel ** 2) + math.sin(sint) ** 2 * (phiVel ** 2)))) / math.sqrt(1 - (r_s / r))

        velocity = np.insert(velocity, 0, timeVel)
        return velocity

    def f(self, s, Y, mass):
        x, v = Y[:4], Y[4:]
        dY_by_ds = np.empty_like(Y)
        dY_by_ds[:4] = v
        (*dvars,) = dt, dr, dθ, dϕ = sp.symbols("dt dr dθ dϕ")
        accel = self.geodesic_equation()
        # because accel is (4,1), it yields one-element lists when iterated
        right = [
            el
            for sublst in accel.subs(
                [(variable, value) for variable, value in zip(self.vars + dvars + [self.M], x.tolist() + v.tolist() + [mass])]
            )
            for el in sublst
        ]
        # print(right)

        dY_by_ds[4:] = right
        # print(Y[4:])
        return dY_by_ds

    def calc_positions(self, to_time, x, v, mass):
        v = self.calc_time_component(x, v)
        res = solve_ivp(
            self.f, t_span=(0, to_time), t_eval=np.linspace(0, to_time, 100), y0=np.array([*x, *v]), args=(mass,),
            method="DOP853"
        )
        return res


class Plot:
    def __init__(self):
        self.sch = Schwarzchild()
        self.initial_pos = np.array([0, 1.5, math.pi / 2, 0])
        self.initial_vel = np.array([0, 0, 0.1])
        for i in range(1, 2):
            new_pos = [0, 1.5 + 0.2 * i, math.pi / 2, 0]
            self.initial_pos = np.vstack((self.initial_pos, new_pos))

        for i in range(1, 2):
            new_vel = [0, 0, 0.1 + 0.1 * i]
            self.initial_vel = np.vstack((self.initial_vel, new_vel))

    def plot(self, mass, to_time):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        for i in range(len(self.initial_pos)):
            position = self.sch.calc_positions(to_time, self.initial_pos[i], self.initial_vel[i], mass)
            proper_times, times, radii, theta, phi = position.t, position.y[0, :], position.y[1, :], position.y[2, :], position.y[3, :]

            X = radii * np.cos(phi)
            Y = radii * np.sin(phi)

            ax.scatter3D(X, Y, times)
        plt.title("X vs Y vs Time")
        plt.show()


animate = Plot()
animate.plot(0.5, 60)
