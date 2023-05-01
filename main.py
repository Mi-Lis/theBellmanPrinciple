# This is a sample Python script.
import sympy as sp
import numpy as np
from scipy.integrate import odeint, simpson
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
k = 1
m = 1
tk = 3
lam = 1
x = [0, 1]


def fs(s, t):
    s1,s2,s3,s4, s5, s6 = s
    return [
        1/4*s3**2,
        1/2*s3*s4,
        s3-s2+s3*s6,
        k / m * s4 - 2 * s5 + s4 * s6,
            1 / 4 * s4 ** 2,
            2 * s6 * k / m - s4 + s6 ** 2]


def f_sys(x, t, s):
    x1, x2 = x
    return [x2, u(x, s, t) - k / m * x2]


def u(x, s, t):
    s3, s4, s6 = s
    x1, x2 = x
    return -1/2*(s4(t) * x1 + 2*s6(t) * x2 + s3(t))


def F3(x):
    x1, x2 = x
    return x1[-1] ** 2 + x2[-1] ** 2


def J(x,s, t):
    return simpson([i ** 2 for i in u(x, s, t)]) + lam * F3(x)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    t = np.linspace(0, tk)

    s456 = odeint(fs, [0,0,0,0, lam, lam], t[::-1])[::-1]

    s3 = s456[:, 2]
    s3 = CubicSpline(t, s3)
    s4 = s456[:, 3]
    s4 = CubicSpline(t, s4)
    s6 = s456[:, 5]
    s6 = CubicSpline(t, s6)
    s = [s3, s4, s6]

    ans = odeint(f_sys, x, t, args=(s,))

    x1 = ans[:, 0]
    x2 = ans[:, 1]

    x = [x1, x2]

    fig5, axs5 = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    axs5[0,0].plot(t, x1, '-c', linewidth=2, label=r"x_1")
    axs5[0,0].legend(loc='best')
    axs5[0,0].xaxis.set_label_coords(1.05, -0.025)
    axs5[0,0].tick_params(labelsize=16)
    axs5[0,0].set_xlabel('$t$', size=16)
    axs5[0,0].set_ylabel(r'$x_1$', rotation=0, size=16)
    axs5[0,0].yaxis.set_label_coords(-.1, .95)
    axs5[0,0].set_title(label="Координата точки")
    axs5[0,0].grid()


    axs5[1,0].plot(t, x2, '-c', linewidth=2, label=r"x_2")
    axs5[1,0].legend(loc='best')
    axs5[1,0].xaxis.set_label_coords(1.05, -0.025)
    axs5[1,0].tick_params(labelsize=16)
    axs5[1,0].set_xlabel('$t$', size=16)
    axs5[1,0].set_ylabel(r'$x_2$', rotation=0, size=16)
    axs5[1,0].yaxis.set_label_coords(-.1, .95)
    axs5[1,0].set_title(label="Скорость точки")
    axs5[1,0].grid()
    plt.show(block=False)


    axs5[0,1].plot(x1, x2, '-c', linewidth=2, label=r"x_1 x_2")
    axs5[0,1].legend(loc='best')
    axs5[0,1].xaxis.set_label_coords(1.05, -0.025)
    axs5[0,1].tick_params(labelsize=16)
    axs5[0,1].set_xlabel(r'$x_1$', size=16)
    axs5[0,1].set_ylabel(r'$x_2$', rotation=0, size=16)
    axs5[0,1].yaxis.set_label_coords(-.1, .95)
    axs5[0,1].set_title(label="Фазовый портрет")
    axs5[0,1].grid()
    plt.show(block=False)

    axs5[1,1].plot(t, u(x, s, t), '-c', linewidth=2, label=r"u")
    axs5[1,1].legend(loc='best')
    axs5[1,1].xaxis.set_label_coords(1.05, -0.025)
    axs5[1,1].tick_params(labelsize=16)
    axs5[1,1].set_xlabel('$t$', size=16)
    axs5[1,1].set_ylabel(r'$x_1$', rotation=0, size=16)
    axs5[1,1].yaxis.set_label_coords(-.1, .95)
    axs5[1,1].set_title(label="Управление движения")
    axs5[1,1].grid()
    plt.show(block=False)

    fig5.savefig('fig/lambda-'+str(lam)+'tk-'+str(tk)+'.png')

    print("J = ", J(x,s, t) ) #100010.00000001
    plt.show()
