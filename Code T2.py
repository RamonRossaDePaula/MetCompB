


# Course: Métodos Computacionais da Física B (FIS01206)


## Note: I decided to write the paper in English as training. I affirm that there may be grammatical errors and I kindly request that they be disregarded.


# Instructions:
'''
Reproduzam os gráficos da página 167 do livro Classical Dynamics of Particles and Systems, (5a ed, aut. Marion, Thornton). 

Use Runge-Kutta 4ª ordem clássico para resolver a equação

    \frac{d^2x}{dt^2}=-c \frac{dx}{dt} - sin (x) + F cos (wt) 

onde \quad x(0) = 1.0, v(0) = 0.0. 

Escreva a equação de 2a ordem como o sistema de equações de 1a ordem.

    \begin{cases} \frac{dx}{dt}=v \\ \frac{dv}{dt}=-c v - sin (x) + F cos (wt)\end{cases} 

com  h=0.01  (x: posição angular, v: velocidade angular). Os outros parâmetros são: c=0.05, w=0.7.
'''

#___________________________________________________________

# Imports

import matplotlib.pyplot as plt
from numpy import pi, sin, cos, linspace
from copy import deepcopy

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