


# Course: Métodos Computacionais da Física B (FIS01206)


## Note: I decided to write the paper in English as training. I affirm that there may be grammatical errors and I kindly request that they be disregarded.


# Instructions:
'''
1) Escolha um problema de valor inicial, de preferência não linear. Não pode ser decaimento/crescimento exponencial!

2) O problema pode ser de primeira ou segunda ordem. O trabalho 2 será um problema de segunda ordem.

3) Resolva o problema usando 3 métodos (ex. Euler explícito, um Runge-Kutta 2 e um Runge-Kutta 4). 
Cuidado: Euler-Cromer e Velocity-Verlet somente podem ser usados se o problema for de segunda ordem.

4) Faça um gráfico mostrando as soluções numéricas e a analítica, se houver.

5) Faça alguma análise do erro: erro de cada método em função do tempo e/ou erro global em função de h. 
Se houver solução analítica, use-a para o cálculo do erro. Se não houver solução analítica, use o resultado 
de Runge-Kutta 4 como a melhor aproximação da solução real e calcule o erro dos outros métodos em relação a este.
'''

#___________________________________________________________

# Imports

import matplotlib.pyplot as plt
from numpy import pi, sin, cos, linspace
from copy import deepcopy

#___________________________________________________________

# IVP:

'''
Resonance case, described by u"+(w**2)u = cos(w0*t)
                          being w = w0 
                          and u'(0) = u(0) = 0

Let's calculate the value of the amplitude at t = 5. For this, we will use w = w0 = π


Source: Whitman College - Course Differential Equations (Math 244 - Spring 2019) - Section 3.8 
Available in: http://people.whitman.edu/~hundledr/courses/M244S19/M244/Sect03_8.pdf
'''

# Definitions:

w = pi

#___________________________________________________________


# Euler-Cromer method:

h = 0.001

def eulercromer(h):
    
    X = [0]
    V = [0]
    t = [0]
    time = 5

    while time > t[-1]:
        vn = V[-1] + (cos(w * t[-1]) - w ** 2 * X[-1]) * h
        xn = X[-1] + vn * h

        V.append(vn)
        X.append(xn)
        t.append(t[-1] + h)

    Xec = deepcopy(X)
    tec = deepcopy(t)

    return tec, Xec

#___________________________________________________________


# Verlet method:


def verlet(h):
  
    X = [0]
    V = [0]
    t = [0]
    time = 5

    vn = V[-1] + (cos(w * t[-1]) - w ** 2 * X[-1]) * h
    xn = X[-1] + vn * h

    V.append(vn)
    X.append(xn)
    t.append(t[-1] + h)

    while time > t[-1]:
        xn = 2 * X[-1] - X[-2] + (cos(w * t[-1]) - w ** 2 * X[-1]) * h ** 2

        V.append(vn)
        X.append(xn)
        t.append(t[-1] + h)

    Xv = deepcopy(X)
    tv = deepcopy(t)

    return tv, Xv

#___________________________________________________________


# Runge-Kutta 4 method:

def rungekutta4(h):
  
    X = [0]
    V = [0]
    t = [0]
    time = 5

    def F(t, x):
        return cos(w * t) - w ** 2 * x

    while t[-1] <= time:
        xm = X[-1]                              
       
        k1_x = V[-1] * h                     
        k1_v = h * F(t[-1], X[-1])
       
        xm = X[-1] + 0.5 * k1_x                  
        vm = V[-1] + 0.5 * k1_v
        k2_x = h * vm
        k2_v = h * F(t[-1] + h / 2, xm)

        xm = X[-1] + 0.5 * k2_x                
        vm = V[-1] + 0.5 * k2_v
        k3_x = h * vm
        k3_v = h * F(t[-1] + h / 2, xm)

        xm = X[-1] + k3_x                      
        vm = V[-1] + k3_v
        k4_x = h * vm
        k4_v = h * F(t[-1] + h, xm)

        x_novo = X[-1] + 1 / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)      
        v_novo = V[-1] + 1 / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

        t.append(t[-1] + h)
        X.append(x_novo)
        V.append(v_novo)

    Xr = deepcopy(X)
    tr = deepcopy(t)

    return tr, Xr

#___________________________________________________________

# Chart settings:

def analytics(h):
    x = linspace(0, 5, num=int(5 / h + 2))
    y = x * sin(pi * x) / (2 * pi)

    return x, y


def __plots__():
    plt.plot(tec, Xec, 'bo', label='Euler-Cromer', color='red')
    plt.plot(tv, Xv, 'bs', markersize=1, label='Velocity-Verlet', color='gold')
    plt.plot(tr, Xr, label='Runge-Kutta 4', color='black')
    plt.plot(t, x, label='Analytical Solution', color='blue')
    plt.legend()
    plt.title('Comparison of methods (h = 0.001)')
    plt.grid()
    plt.show()


tec, Xec = eulercromer(h)
tv, Xv = verlet(h)
tr, Xr = rungekutta4(h)
t, x = analytics(h)

#___________________________________________________________

# Error:

def __error__():
    ERROR_EULERCROMER = []
    ERROR_VERLET = []
    ERROR_RUNGEKUTTA = []

    for i in range(5001):
        error_ec = abs(Xec[i] - x[i])
        error_v = abs(Xv[i] - x[i])
        error_r = abs(Xr[i] - x[i])

        ERROR_EULERCROMER.append(error_ec)
        ERROR_VERLET.append(error_v)
        ERROR_RUNGEKUTTA.append(error_r)

    plt.plot(tec,  ERROR_EULERCROMER, 'o', markersize=5, color='red', label='Error: Euler-Cromer')
    plt.plot(tv, ERROR_VERLET, 'bs', markersize=3, color='gold', label='Error: Verlet')
    plt.plot(tr, ERROR_RUNGEKUTTA, 'x', markersize=1, color='black', label='Error: Runge-Kutta 4')

    plt.yscale('log')
    plt.title('Error analysis (h = 0.001)')
    plt.legend()
    plt.grid()
    plt.show()


__plots__()
__error__()