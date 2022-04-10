


# Course: Métodos Computacionais da Física B (FIS01206)

# Test 1

#___________________________________________________________

# Question 1

# Imports
from numpy import pi
import matplotlib.pyplot as plt

#___________________________________________________________

# Data

# dx/dt = -0.6 * pi * (r**2) * ((2 * g * x)^0.5)/A(x)
r = 0.1     # pé
g = 32.1    # pés/s²
# k1 = f(tn, xn)
# k2 = f(tn + h, xn + k1 * h)
# xn+1 = xn + (h/2) * (k1 + k2)

#___________________________________________________________


# Functions

def A(x):
    return pi * x ** 2


def f(xn):
    return -0.6 * pi * (r ** 2) * ((2 * g * xn) ** 0.5) / A(xn)


def flowrate(letter='b', tf=600, h=0.1):
    X = [8]                 # feet
    T = [0]
    VOLUME = [536.16]       # feet³

    if letter == 'b':

        while tf > T[-1]:
            k1 = f(X[-1])
            k2 = f(X[-1] + k1 * h)

            xn = X[-1] + (h / 2) * (k1 + k2)

            volume = (pi / 3) * xn ** 3

            VOLUME.append(volume)
            X.append(xn)
            T.append(T[-1] + h)

    if letter == 'c':
        condition = True

        while condition:

            k1 = f(X[-1])
            k2 = f(X[-1] + k1 * h)

            xn = X[-1] + (h / 2) * (k1 + k2)

            volume = (pi / 3) * xn ** 3

            VOLUME.append(volume)
            X.append(xn)
            T.append(T[-1] + h)

            if X[-1] < 0.1:
                condition = False

    return X, T, VOLUME


def graphics():

    X, T, V = flowrate('b')

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs[0, 0].plot(T, X, color='crimson')
    axs[0, 0].set_title('Radius varying in 10 min.')
    axs[0, 0].set_xlabel('Time (s)')
    axs[0, 0].set_ylabel('Radius (feet)')

    axs[0, 1].plot(T, V)
    axs[0, 1].set_title(f'Volume at the instant of 10 min.: {round(V[-1], 2)}feet³')
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Volume (feet³)')

    X, T, V = flowrate('c')

    axs[1, 0].plot(T, X, color='crimson')
    axs[1, 0].set_title(f'Radius varying in {round(T[-1]/60, 2)} minutes.')
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_ylabel('Radius (feet)')

    axs[1, 1].plot(T, V)
    axs[1, 1].set_title(f'Volume at the instant of {round(T[-1]/60, 2)} minutes: {round(V[-1], 2)}feet³')
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_ylabel('Volume (feet³)')

    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.2,
                        hspace=0.4)

    plt.show()


graphics()