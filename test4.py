import sympy as sp, numpy as np
import matplotlib.animation as animation
import math

t, r, θ, ϕ, M = sp.symbols("t r θ ϕ M")
vars = [t, r, θ, ϕ]
g_co = sp.Matrix(
    [
        [-(1 - 2 * M / r), 0, 0, 0],
        [0, 1 / (1 - 2 * M / r), 0, 0],
        [0, 0, r ** 2, 0],
        [0, 0, 0, r ** 2 * sp.sin(θ) ** 2],
    ]
)

g_contra = g_co.inv()

dg_co_by = sp.Array([[[sp.diff(g_co[i, j], vars[k]) for k in range(4)] for j in range(4)] for i in range(4)])


def basic_sum(it):
    it = iter(it)
    cur = next(it)
    for el in it:
        cur += el
    return cur


Γ = sp.Array(
    [
        [
            [
                basic_sum(
                    g_contra[i, m] * (dg_co_by[m, k, l] + dg_co_by[m, l, k] - dg_co_by[k, l, m]) for m in range(4)
                )
                / 2
                for l in range(4)
            ]
            for k in range(4)
        ]
        for i in range(4)
    ]
)

# print(Γ)

s = sp.symbols("s")
(*dvars,) = dt, dr, dθ, dϕ = sp.symbols("dt dr dθ dϕ")  # "speeds" (derivatives by ds) of each coord
accel = sp.Array(
    [
        -basic_sum(Γ[mu, alpha, beta] * (dvars[alpha] * dvars[beta]) for alpha in range(4) for beta in range(4))
        for mu in range(4)
    ],
    shape=(4, 1),
)

# print(accel)

from scipy.integrate import solve_ivp


def f(s, Y, mass):
    x, v = Y[:4], Y[4:]
    dY_by_ds = np.empty_like(Y)
    dY_by_ds[:4] = v
    # because accel is (4,1), it yields one-element lists when iterated
    right = [
        el
        for sublst in accel.subs(
            [(variable, value) for variable, value in zip(vars + dvars + [M], x.tolist() + v.tolist() + [mass])]
        )
        for el in sublst
    ]
    # print(right)
    dY_by_ds[4:] = right
    return dY_by_ds


mass = 1
r0 = 5 * mass
x = np.array([0, r0, np.pi / 2, 0])
v = np.array([1 / np.sqrt(1 - 2 * mass / r0), 0, 0, 0])
to_time = 11


res = solve_ivp(
    f, t_span=(0, to_time), t_eval=np.linspace(0, to_time, 100), y0=np.array([*x, *v]), args=(mass,), method="DOP853"
)

proper_times, times, radii = res.t, res.y[0, :], res.y[1, :]

import matplotlib.pyplot as plt

fig, [ax1, ax2] = plt.subplots(1, 2)
ax1.plot(proper_times, times, label="t(s)")
ax2.plot(times, radii, label="r(s)")
for ax in [ax1, ax2]:
    ax.legend()
    ax.set_xlabel("s (proper time)")
plt.show()

