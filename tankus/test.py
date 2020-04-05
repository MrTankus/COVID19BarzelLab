import numpy as np
from scipy.integrate import odeint


def model(z, t, alpha, beta, gama, population):
    s = z[0]
    e = z[1]
    i = z[2]
    dsdt = -beta * s * (e / population)
    dedt = beta * e * (s / population) - gama * e
    didt = gama * e - alpha * i
    drdt = alpha * i
    return [dsdt, dedt, didt, drdt]


days = 300
total_hours = 24 * days

z0 = [0, 0, 0]
t = np.linspace(0, days, total_hours)
alpha = 0.004
beta = 0.001
gama = 0.00083
population = 9000000

ventilators = 3000
need_ventilators = 0.01

z = odeint(model, z0, t, args=())

s = z[:, 0]
e = z[:, 1]
i = z[:, 2]
r = z[:, 3]